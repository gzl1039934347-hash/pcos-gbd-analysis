# Nordpred预测PCOS发病人数和年龄标化发病率

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c(
  "nordpred", "dplyr", "tidyr", "readr", "ggplot2", "forecast", "patchwork"
))

NORDPRED_BOOTSTRAP_REPS <- 500L
OBS_PERIOD_STARTS <- seq(1994L, 2019L, by = 5L)
FORECAST_PERIOD_STARTS <- seq(2024L, 2034L, by = 5L)
OBS_PERIOD_LABELS <- paste0(OBS_PERIOD_STARTS, "-", OBS_PERIOD_STARTS + 4L)
FORECAST_PERIOD_LABELS <- paste0(FORECAST_PERIOD_STARTS, "-", FORECAST_PERIOD_STARTS + 4L)
FORECAST_PERIOD_MID <- FORECAST_PERIOD_STARTS + 2L

pcos <- readRDS(PCOS_CLEAN_RDS)
population_historical <- readRDS(POPULATION_RDS)
population_future <- readRDS(WPP_RDS)
weights <- readr::read_csv(AGE_WEIGHT_CSV, show_col_types = FALSE) |>
  dplyr::arrange(match(age_name, AGE_GROUPS))

incidence <- pcos |>
  dplyr::filter(
    measure == "Incidence", metric_name == "Number",
    age_name %in% AGE_GROUPS
  ) |>
  dplyr::select(location_id, location_name, year, age_name, val, lower, upper)

period_label <- function(year, starts) {
  hit <- starts[findInterval(year, starts)]
  ifelse(year >= min(starts) & year <= max(starts) + 4L,
         paste0(hit, "-", hit + 4L), NA_character_)
}

make_group_annual <- function(ids, group_id, group_name) {
  cases <- incidence |>
    dplyr::filter(location_id %in% ids) |>
    dplyr::group_by(year, age_name) |>
    dplyr::summarise(
      val = sum(val), lower = sum(lower), upper = sum(upper), .groups = "drop"
    ) |>
    dplyr::mutate(location_id = group_id, location_name = group_name)
  population <- dplyr::bind_rows(population_historical, population_future) |>
    dplyr::filter(location_id %in% ids) |>
    dplyr::group_by(year, age_name) |>
    dplyr::summarise(population = sum(population), .groups = "drop") |>
    dplyr::mutate(location_id = group_id, location_name = group_name)
  list(cases = cases, population = population)
}

country_groups <- lapply(seq_len(nrow(COUNTRIES)), function(i) {
  make_group_annual(
    COUNTRIES$location_id[[i]], COUNTRIES$location_id[[i]], COUNTRIES$location_name[[i]]
  )
})
names(country_groups) <- as.character(COUNTRIES$location_id)
country_groups[["-1"]] <- make_group_annual(
  setdiff(COUNTRIES$location_id, 6L), -1L, "14 neighbouring countries combined"
)

aggregate_periods <- function(group, observed_starts, future_starts) {
  observed <- group$cases |>
    dplyr::filter(year >= min(observed_starts), year <= max(observed_starts) + 4L) |>
    dplyr::mutate(period = period_label(year, observed_starts)) |>
    dplyr::group_by(age_name, period) |>
    dplyr::summarise(
      val = sum(val), lower = sum(lower), upper = sum(upper), .groups = "drop"
    )
  pop <- group$population |>
    dplyr::filter(
      year >= min(observed_starts),
      year <= max(c(observed_starts, future_starts)) + 4L
    ) |>
    dplyr::mutate(
      period = dplyr::if_else(
        year <= max(observed_starts) + 4L,
        period_label(year, observed_starts),
        period_label(year, future_starts)
      )
    ) |>
    dplyr::filter(!is.na(period)) |>
    dplyr::group_by(age_name, period) |>
    dplyr::summarise(population = sum(population), .groups = "drop")
  list(cases = observed, population = pop)
}

to_matrix <- function(data, value, periods) {
  data |>
    dplyr::mutate(
      age_name = factor(age_name, levels = AGE_GROUPS),
      period = factor(period, levels = periods)
    ) |>
    dplyr::select(age_name, period, value = dplyr::all_of(value)) |>
    tidyr::complete(age_name, period) |>
    tidyr::pivot_wider(names_from = period, values_from = value) |>
    dplyr::arrange(age_name) |>
    tibble::column_to_rownames("age_name") |>
    as.data.frame()
}

# 官方nordpred R实现将年龄组数量硬编码为18；以下兼容函数仅将该常数
# 替换为输入矩阵的实际年龄组数，其余模型、链接函数和趋势截断逻辑不变。
nordpred_flexible_estimate <- function(cases, pyr, noperiod, startestage,
                                       linkfunc = "power5") {
  dnoagegr <- nrow(cases)
  if (nrow(pyr) != dnoagegr) stop("cases与pyr年龄组数量不一致。")
  if (ncol(cases) > ncol(pyr)) stop("pyr必须覆盖全部观察期。")
  if (ncol(pyr) == ncol(cases)) stop("pyr必须包含预测期人口。")
  if ((ncol(pyr) - ncol(cases)) > 5) stop("Nordpred最多预测5个时期。")
  if (noperiod < 3 || noperiod > ncol(cases)) stop("noperiod必须为3至观察期数。")
  dnoperiods <- ncol(cases)
  ageno <- rep(seq_len(dnoagegr), dnoperiods)
  periodno <- sort(rep(seq_len(dnoperiods), dnoagegr))
  cohort <- max(ageno) - ageno + periodno
  y <- c(as.matrix(pyr[, seq_len(dnoperiods), drop = FALSE]))
  apcdata <- data.frame(
    Cases = c(as.matrix(cases)), Age = ageno, Cohort = cohort,
    Period = periodno, y = y
  ) |>
    dplyr::filter(Age >= startestage, Period > dnoperiods - noperiod)

  power5link <- stats::poisson()
  power5link$link <- "0.2 root link Poisson family"
  power5link$linkfun <- function(mu) (mu / apcdata$y)^0.2
  power5link$linkinv <- function(eta) pmax(.Machine$double.eps, apcdata$y * eta^5)
  power5link$mu.eta <- function(eta) pmax(.Machine$double.eps, 5 * apcdata$y * eta^4)
  old_contrasts <- getOption("contrasts")
  on.exit(options(contrasts = old_contrasts), add = TRUE)
  options(contrasts = c("contr.treatment", "contr.poly"))

  if (linkfunc == "power5") {
    model <- stats::glm(
      Cases ~ as.factor(Age) + Period + as.factor(Period) +
        as.factor(Cohort) - 1,
      family = power5link, data = apcdata
    )
  } else if (linkfunc == "poisson") {
    model <- stats::glm(
      Cases ~ as.factor(Age) + Period + as.factor(Period) +
        as.factor(Cohort) + offset(log(y)) - 1,
      family = stats::poisson(), data = apcdata
    )
  } else {
    stop("未知linkfunc。")
  }
  pvalue <- 1 - stats::pchisq(model$deviance, model$df.residual)
  mod1 <- stats::glm(
    Cases ~ as.factor(Age) + Period + as.factor(Cohort) + offset(log(y)) - 1,
    family = stats::poisson(), data = apcdata
  )
  mod2 <- stats::glm(
    Cases ~ as.factor(Age) + Period + I(Period^2) + as.factor(Cohort) +
      offset(log(y)) - 1,
    family = stats::poisson(), data = apcdata
  )
  pdiff <- stats::anova(mod1, mod2, test = "Chisq")[["Pr(>Chi)"]][2]
  if (is.null(pdiff)) pdiff <- stats::anova(mod1, mod2, test = "Chisq")[["P(>|Chi|)"]][2]
  result <- list(
    glm = model, cases = cases, pyr = pyr, noperiod = noperiod,
    gofpvalue = pvalue, startestage = startestage,
    suggestionrecent = isTRUE(pdiff < 0.05), pvaluerecent = pdiff,
    linkfunc = linkfunc
  )
  class(result) <- "nordpred.estimate"
  result
}

nordpred_flexible_prediction <- function(estimate, startuseage, recent, cuttrend) {
  cases <- as.matrix(estimate$cases)
  pyr <- as.matrix(estimate$pyr)
  dnoagegr <- nrow(cases)
  noperiod <- estimate$noperiod
  nototper <- ncol(pyr)
  noobsper <- ncol(cases)
  nonewpred <- nototper - noobsper
  cuttrend <- cuttrend[seq_len(nonewpred)]
  years <- colnames(pyr)
  if (is.null(years)) years <- paste0("Period", seq_len(nototper))
  datatable <- matrix(NA_real_, dnoagegr, nototper)
  datatable[, seq_len(noobsper)] <- as.matrix(cases)
  rownames(datatable) <- rownames(cases)
  colnames(datatable) <- years
  if (startuseage > 1) {
    for (age in seq_len(startuseage - 1L)) {
      obsinc <- cases[age, (noobsper - 1L):noobsper] /
        pyr[age, (noobsper - 1L):noobsper]
      obsinc[is.na(obsinc)] <- 0
      datatable[age, (noobsper + 1L):nototper] <-
        mean(obsinc) * pyr[age, (noobsper + 1L):nototper]
    }
  }
  for (age in startuseage:dnoagegr) {
    startestage <- estimate$startestage
    coefficients <- estimate$glm$coefficients
    coh <- (dnoagegr - startestage) - (age - startestage) +
      (noperiod + seq_len(nonewpred))
    noages <- dnoagegr - startestage + 1L
    driftmp <- cumsum(1 - cuttrend)
    cohfind <- noages + (noperiod - 1L) + 1L + (coh - 1L)
    maxcoh <- dnoagegr - startuseage + noperiod
    agepar <- as.numeric(coefficients[age - startestage + 1L])
    driftfind <- pmatch("Period", names(coefficients))
    driftpar <- as.numeric(coefficients[driftfind])
    cohpar <- vapply(coh, function(current_cohort) {
      if (current_cohort < maxcoh) {
        as.numeric(coefficients[cohfind[which(coh == current_cohort)[1]]])
      } else {
        value <- as.numeric(coefficients[length(coefficients) -
                                           (startuseage - startestage)])
        ifelse(is.na(value), 0, value)
      }
    }, numeric(1))
    if (recent) {
      lpfind <- driftfind + noperiod - 2L
      driftrecent <- driftpar - as.numeric(coefficients[lpfind])
    }
    if (estimate$linkfunc == "power5") {
      rate <- if (recent) {
        (agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)^5
      } else {
        (agepar + driftpar * (noobsper + driftmp) + cohpar)^5
      }
    } else {
      rate <- if (recent) {
        exp(agepar + driftpar * noobsper + driftrecent * driftmp + cohpar)
      } else {
        exp(agepar + driftpar * (noobsper + driftmp) + cohpar)
      }
    }
    datatable[age, (noobsper + 1L):nototper] <-
      rate * pyr[age, (noobsper + 1L):nototper]
  }
  result <- list(
    predictions = as.data.frame(datatable), pyr = as.data.frame(pyr),
    cases = as.data.frame(cases), linkfunc = estimate$linkfunc,
    nopred = nonewpred, noperiod = noperiod, gofpvalue = estimate$gofpvalue,
    recent = recent, pvaluerecent = estimate$pvaluerecent,
    cuttrend = cuttrend, startuseage = startuseage,
    startestage = estimate$startestage, glm = estimate$glm
  )
  class(result) <- "nordpred"
  result
}

nordpred_flexible <- function(cases, pyr, startestage, startuseage,
                              noperiods, recent = NULL, cuttrend,
                              linkfunc = "power5") {
  candidates <- sort(noperiods)
  while (length(candidates) > 1) {
    widest <- max(candidates)
    model <- nordpred_flexible_estimate(cases, pyr, widest, startestage, linkfunc)$glm
    pvalue <- 1 - stats::pchisq(model$deviance, model$df.residual)
    candidates <- if (pvalue < 0.01) candidates[-length(candidates)] else widest
  }
  noperiod <- candidates[[1]]
  estimate <- nordpred_flexible_estimate(cases, pyr, noperiod, startestage, linkfunc)
  if (is.null(recent)) recent <- estimate$suggestionrecent
  nordpred_flexible_prediction(estimate, startuseage, recent, cuttrend)
}

fit_nordpred <- function(cases_matrix, population_matrix, future_labels,
                         linkfunc = "power5") {
  nordpred_flexible(
    cases = as.data.frame(round(cases_matrix)),
    pyr = as.data.frame(population_matrix),
    startestage = 1, startuseage = 1,
    noperiods = 4:ncol(cases_matrix), recent = NULL,
    cuttrend = c(0, 0.25, 0.50)[seq_along(future_labels)],
    linkfunc = linkfunc
  )
}

extract_nordpred <- function(fit, location_id, location_name, future_labels) {
  counts <- as.matrix(nordpred::nordpred.getpred(
    fit, incidence = FALSE, excludeobs = TRUE, byage = TRUE
  ))
  rates <- as.matrix(nordpred::nordpred.getpred(
    fit, incidence = TRUE, excludeobs = TRUE, byage = TRUE
  ))
  asr <- as.numeric(nordpred::nordpred.getpred(
    fit, incidence = TRUE, standpop = weights$weight,
    excludeobs = TRUE, byage = FALSE
  ))
  colnames(counts) <- future_labels
  colnames(rates) <- future_labels
  age <- expand.grid(
    age_name = rownames(counts), period = colnames(counts),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      location_id = location_id,
      location_name = location_name,
      number = as.vector(counts),
      rate = as.vector(rates)
    ) |>
    dplyr::select(location_id, location_name, period, age_name, number, rate)
  total <- data.frame(
    location_id = location_id, location_name = location_name,
    period = future_labels, period_mid = as.integer(substr(future_labels, 1, 4)) + 2L,
    number = colSums(counts), asr = asr
  )
  list(age = age, total = total)
}

sample_cases_matrix <- function(period_data, observed_labels) {
  sampled <- period_data |>
    dplyr::mutate(
      log_se = dplyr::if_else(
        lower > 0 & upper > lower,
        (log(upper) - log(lower)) / (2 * 1.96), 0
      ),
      sampled = dplyr::if_else(
        log_se > 0,
        stats::rlnorm(dplyr::n(), meanlog = log(pmax(val, 1e-8)), sdlog = log_se),
        val
      )
    )
  as.matrix(to_matrix(sampled, "sampled", observed_labels))
}

bootstrap_forecast <- function(period_data, population_matrix, observed_labels,
                               future_labels, linkfunc = "power5",
                               reps = NORDPRED_BOOTSTRAP_REPS,
                               seed = RANDOM_SEED) {
  set.seed(seed)
  cases_matrix <- as.matrix(to_matrix(period_data, "val", observed_labels))
  base_fit <- fit_nordpred(
    cases_matrix, population_matrix, future_labels, linkfunc = linkfunc
  )
  covariance <- stats::vcov(base_fit$glm)
  coefficient_names <- intersect(rownames(covariance), names(stats::coef(base_fit$glm)))
  covariance <- covariance[coefficient_names, coefficient_names, drop = FALSE]
  coefficient_mean <- stats::coef(base_fit$glm)[coefficient_names]
  estimable <- is.finite(coefficient_mean) & is.finite(diag(covariance))
  coefficient_names <- coefficient_names[estimable]
  coefficient_mean <- coefficient_mean[estimable]
  covariance <- covariance[estimable, estimable, drop = FALSE]
  complete_covariance <- apply(covariance, 1, function(x) all(is.finite(x)))
  coefficient_names <- coefficient_names[complete_covariance]
  coefficient_mean <- coefficient_mean[complete_covariance]
  covariance <- covariance[complete_covariance, complete_covariance, drop = FALSE]
  eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  eig$values <- pmax(eig$values, max(eig$values) * 1e-10)
  covariance <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  dimnames(covariance) <- list(coefficient_names, coefficient_names)
  coefficient_draws <- MASS::mvrnorm(
    reps,
    mu = coefficient_mean,
    Sigma = covariance
  )
  if (reps == 1L) coefficient_draws <- matrix(coefficient_draws, nrow = 1L)
  colnames(coefficient_draws) <- coefficient_names
  draws <- lapply(seq_len(reps), function(i) {
    tryCatch({
      estimate <- list(
        glm = base_fit$glm, cases = base_fit$cases, pyr = base_fit$pyr,
        noperiod = base_fit$noperiod, gofpvalue = base_fit$gofpvalue,
        startestage = base_fit$startestage,
        suggestionrecent = base_fit$recent,
        pvaluerecent = base_fit$pvaluerecent,
        linkfunc = base_fit$linkfunc
      )
      estimate$glm$coefficients[coefficient_names] <- coefficient_draws[i, ]
      class(estimate) <- "nordpred.estimate"
      fit <- nordpred_flexible_prediction(
        estimate, base_fit$startuseage, base_fit$recent, base_fit$cuttrend
      )
      data.frame(
        period = future_labels,
        number = colSums(as.matrix(nordpred::nordpred.getpred(
          fit, incidence = FALSE, excludeobs = TRUE, byage = TRUE
        ))),
        asr = as.numeric(nordpred::nordpred.getpred(
          fit, incidence = TRUE, standpop = weights$weight,
          excludeobs = TRUE, byage = FALSE
        ))
      )
    }, error = function(e) NULL)
  }) |>
    dplyr::bind_rows(.id = "draw")
  assert_true(length(unique(draws$draw)) >= ceiling(0.8 * reps),
              "Nordpred bootstrap成功率低于80%。")
  draws |>
    dplyr::group_by(period) |>
    dplyr::summarise(
      number_lower = stats::quantile(number, 0.025, na.rm = TRUE),
      number_upper = stats::quantile(number, 0.975, na.rm = TRUE),
      asr_lower = stats::quantile(asr, 0.025, na.rm = TRUE),
      asr_upper = stats::quantile(asr, 0.975, na.rm = TRUE),
      successful_draws = dplyr::n_distinct(draw), .groups = "drop"
    )
}

forecast_ids <- c(COUNTRIES$location_id, -1L)
model_specifications <- data.frame(
  model = c("Power5 APC", "Poisson APC"),
  linkfunc = c("power5", "poisson"),
  stringsAsFactors = FALSE
)

annual_incidence <- function(group, end_year = YEAR_END) {
  group$cases |>
    dplyr::filter(year <= end_year) |>
    dplyr::inner_join(
      group$population |>
        dplyr::filter(year <= end_year) |>
        dplyr::select(year, age_name, population),
      by = c("year", "age_name")
    ) |>
    dplyr::inner_join(dplyr::select(weights, age_name, weight), by = "age_name") |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      number = sum(val),
      asr = sum(val / population * RATE_SCALE * weight),
      .groups = "drop"
    ) |>
    dplyr::arrange(year)
}

fit_ets_outcome <- function(annual, outcome, train_end, forecast_end) {
  training <- annual |>
    dplyr::filter(year <= train_end) |>
    dplyr::arrange(year)
  assert_true(max(training$year) == train_end, "ETS训练截止年份不完整。")
  fit <- forecast::ets(stats::ts(
    training[[outcome]], start = min(training$year), frequency = 1
  ))
  fc <- forecast::forecast(fit, h = forecast_end - train_end, level = 95)
  data.frame(
    year = (train_end + 1L):forecast_end,
    estimate = pmax(as.numeric(fc$mean), 0),
    lower = pmax(as.numeric(fc$lower[, 1]), 0),
    upper = pmax(as.numeric(fc$upper[, 1]), 0)
  )
}

aggregate_ets_periods <- function(number_fc, asr_fc, period_starts) {
  labels <- paste0(period_starts, "-", period_starts + 4L)
  dplyr::left_join(
    number_fc |>
      dplyr::mutate(period = period_label(year, period_starts)) |>
      dplyr::filter(period %in% labels) |>
      dplyr::group_by(period) |>
      dplyr::summarise(
        number = sum(estimate), number_lower = sum(lower),
        number_upper = sum(upper), .groups = "drop"
      ),
    asr_fc |>
      dplyr::mutate(period = period_label(year, period_starts)) |>
      dplyr::filter(period %in% labels) |>
      dplyr::group_by(period) |>
      dplyr::summarise(
        asr = mean(estimate), asr_lower = mean(lower),
        asr_upper = mean(upper), .groups = "drop"
      ),
    by = "period"
  ) |>
    dplyr::mutate(period_mid = as.integer(substr(period, 1, 4)) + 2L)
}

observed_periods <- function(group, period_starts) {
  labels <- paste0(period_starts, "-", period_starts + 4L)
  annual_incidence(group) |>
    dplyr::mutate(period = period_label(year, period_starts)) |>
    dplyr::filter(period %in% labels) |>
    dplyr::group_by(period) |>
    dplyr::summarise(
      actual_number = sum(number), actual_asr = mean(asr), .groups = "drop"
    )
}

run_nordpred <- function(group, location_id, location_name, observed_starts,
                         future_starts, model, linkfunc, seed) {
  observed_labels <- paste0(observed_starts, "-", observed_starts + 4L)
  future_labels <- paste0(future_starts, "-", future_starts + 4L)
  pd <- aggregate_periods(group, observed_starts, future_starts)
  cases_matrix <- as.matrix(to_matrix(pd$cases, "val", observed_labels))
  population_matrix <- as.matrix(to_matrix(
    pd$population, "population", c(observed_labels, future_labels)
  ))
  assert_true(!anyNA(cases_matrix) && !anyNA(population_matrix),
              paste0(location_name, "的Nordpred矩阵存在缺失。"))
  fit <- fit_nordpred(
    cases_matrix, population_matrix, future_labels, linkfunc = linkfunc
  )
  point <- extract_nordpred(fit, location_id, location_name, future_labels)
  interval <- bootstrap_forecast(
    pd$cases, population_matrix, observed_labels, future_labels,
    linkfunc = linkfunc, reps = NORDPRED_BOOTSTRAP_REPS, seed = seed
  )
  point$total <- point$total |>
    dplyr::left_join(interval, by = "period") |>
    dplyr::mutate(model = model)
  point$age <- point$age |>
    dplyr::mutate(model = model)
  point
}

top3_ids <- incidence |>
  dplyr::filter(year == YEAR_END, location_id != 6L) |>
  dplyr::group_by(location_id, location_name) |>
  dplyr::summarise(cases = sum(val), .groups = "drop") |>
  dplyr::slice_max(cases, n = 3, with_ties = FALSE) |>
  dplyr::pull(location_id)
selected_ids <- unique(c(6L, -1L, top3_ids))

future_nordpred <- list()
future_age <- list()
future_ets_annual <- list()
counter <- 1L
for (id in forecast_ids) {
  group <- country_groups[[as.character(id)]]
  name <- unique(group$population$location_name)[[1]]
  for (model_index in seq_len(nrow(model_specifications))) {
    specification <- model_specifications[model_index, ]
    message("未来预测: ", name, " / ", specification$model)
    result <- run_nordpred(
      group, id, name, OBS_PERIOD_STARTS, FORECAST_PERIOD_STARTS,
      specification$model, specification$linkfunc,
      RANDOM_SEED + abs(id) * 10L + model_index
    )
    future_nordpred[[counter]] <- result$total
    future_age[[counter]] <- result$age
    counter <- counter + 1L
  }
  annual <- annual_incidence(group)
  number_fc <- fit_ets_outcome(annual, "number", YEAR_END, FORECAST_END)
  asr_fc <- fit_ets_outcome(annual, "asr", YEAR_END, FORECAST_END)
  future_ets_annual[[as.character(id)]] <- dplyr::left_join(
    dplyr::rename(number_fc, number = estimate, number_lower = lower,
                  number_upper = upper),
    dplyr::rename(asr_fc, asr = estimate, asr_lower = lower, asr_upper = upper),
    by = "year"
  ) |>
    dplyr::mutate(location_id = id, location_name = name, model = "ETS")
}

future_ets_annual <- dplyr::bind_rows(future_ets_annual)
future_ets <- future_ets_annual |>
  dplyr::group_by(location_id, location_name) |>
  dplyr::group_modify(~ aggregate_ets_periods(
    dplyr::select(.x, year, estimate = number, lower = number_lower,
                  upper = number_upper),
    dplyr::select(.x, year, estimate = asr, lower = asr_lower,
                  upper = asr_upper),
    FORECAST_PERIOD_STARTS
  )) |>
  dplyr::ungroup() |>
  dplyr::mutate(model = "ETS")

future_forecast <- dplyr::bind_rows(
  dplyr::bind_rows(future_nordpred), future_ets
) |>
  dplyr::select(
    location_id, location_name, model, period, period_mid,
    number, number_lower, number_upper, asr, asr_lower, asr_upper,
    dplyr::any_of("successful_draws")
  ) |>
  dplyr::arrange(location_id, model, period_mid)
assert_true(nrow(future_forecast) == 16L * 3L * 3L,
            "未来预测汇总记录数不等于144。")
assert_true(!anyNA(dplyr::select(
  future_forecast, number, number_lower, number_upper, asr, asr_lower, asr_upper
)), "未来预测存在缺失值。")
assert_true(all(future_forecast$number_lower <= future_forecast$number_upper) &&
              all(future_forecast$asr_lower <= future_forecast$asr_upper),
            "未来预测区间顺序错误。")

save_table(future_forecast, "Supplement_forecast_future_all_models.csv")
save_table(future_ets_annual, "Supplement_ETS_annual_2024_2040.csv")
save_table(
  dplyr::bind_rows(future_age),
  "Supplement_Nordpred_age_specific_forecast.csv"
)

validation_folds <- list(
  list(name = "origin2013_h5", train_starts = seq(1994L, 2009L, by = 5L),
       prediction_starts = 2014L, origin = 2013L),
  list(name = "origin2013_h10", train_starts = seq(1994L, 2009L, by = 5L),
       prediction_starts = c(2014L, 2019L), origin = 2013L),
  list(name = "origin2018_h5", train_starts = seq(1994L, 2014L, by = 5L),
       prediction_starts = 2019L, origin = 2018L)
)

validation_predictions <- list()
counter <- 1L
for (id in forecast_ids) {
  group <- country_groups[[as.character(id)]]
  name <- unique(group$population$location_name)[[1]]
  annual <- annual_incidence(group)
  for (fold_index in seq_along(validation_folds)) {
    validation_spec <- validation_folds[[fold_index]]
    fold_name <- validation_spec$name
    train_starts <- validation_spec$train_starts
    prediction_starts <- validation_spec$prediction_starts
    origin_year <- validation_spec$origin
    observed <- observed_periods(group, prediction_starts)
    for (model_index in seq_len(nrow(model_specifications))) {
      specification <- model_specifications[model_index, ]
      message("滚动验证: ", name, " / ", specification$model, " / ", fold_name)
      result <- run_nordpred(
        group, id, name, train_starts, prediction_starts,
        specification$model, specification$linkfunc,
        RANDOM_SEED + abs(id) * 100L + fold_index * 10L + model_index
      )$total
      validation_predictions[[counter]] <- result |>
        dplyr::left_join(observed, by = "period") |>
        dplyr::mutate(
          fold = fold_name, origin_year = origin_year,
          training_end = max(train_starts) + 4L
        )
      counter <- counter + 1L
    }
    forecast_end <- max(prediction_starts) + 4L
    ets_number <- fit_ets_outcome(annual, "number", origin_year, forecast_end)
    ets_asr <- fit_ets_outcome(annual, "asr", origin_year, forecast_end)
    ets_result <- aggregate_ets_periods(
      ets_number, ets_asr, prediction_starts
    ) |>
      dplyr::mutate(
        location_id = id, location_name = name, model = "ETS",
        fold = fold_name, origin_year = origin_year,
        training_end = origin_year
      ) |>
      dplyr::left_join(observed, by = "period")
    validation_predictions[[counter]] <- ets_result
    counter <- counter + 1L
  }
}

validation_predictions <- dplyr::bind_rows(validation_predictions) |>
  dplyr::select(
    location_id, location_name, model, fold, origin_year, training_end,
    period, period_mid, actual_number, number, number_lower, number_upper,
    actual_asr, asr, asr_lower, asr_upper, dplyr::any_of("successful_draws")
  ) |>
  dplyr::arrange(location_id, model, fold, period_mid)
assert_true(all(validation_predictions$training_end <= validation_predictions$origin_year),
            "验证训练集包含了预测期数据。")
assert_true(!anyNA(dplyr::select(
  validation_predictions, actual_number, number, number_lower, number_upper,
  actual_asr, asr, asr_lower, asr_upper
)), "验证预测存在缺失值。")
save_table(validation_predictions, "Supplement_forecast_validation_predictions.csv")

validation_long <- dplyr::bind_rows(
  validation_predictions |>
    dplyr::transmute(
      location_id, location_name, model, fold, period_mid, outcome = "Incidence number",
      actual = actual_number, predicted = number,
      lower = number_lower, upper = number_upper
    ),
  validation_predictions |>
    dplyr::transmute(
      location_id, location_name, model, fold, period_mid,
      outcome = "Age-standardized incidence rate",
      actual = actual_asr, predicted = asr,
      lower = asr_lower, upper = asr_upper
    )
)

trend_direction <- function(actual, predicted, fold, period_mid) {
  keep <- fold == "origin2013_h10"
  if (sum(keep) < 2L) return(NA)
  ordering <- order(period_mid[keep])
  sign(diff(actual[keep][ordering])) == sign(diff(predicted[keep][ordering]))
}

model_metrics <- validation_long |>
  dplyr::group_by(location_id, location_name, model, outcome) |>
  dplyr::summarise(
    MAPE = mean(abs(predicted - actual) / pmax(abs(actual), 1e-12)) * 100,
    sMAPE = mean(2 * abs(predicted - actual) /
                   pmax(abs(actual) + abs(predicted), 1e-12)) * 100,
    MAE = mean(abs(predicted - actual)),
    RMSE = sqrt(mean((predicted - actual)^2)),
    coverage_95 = mean(actual >= lower & actual <= upper) * 100,
    mean_interval_width = mean(upper - lower),
    trend_consistent = trend_direction(actual, predicted, fold, period_mid),
    validation_records = dplyr::n(), .groups = "drop"
  )
assert_true(nrow(model_metrics) == 16L * 3L * 2L,
            "验证指标未覆盖16个地区、3个模型和2个结局。")

model_ranking <- model_metrics |>
  dplyr::group_by(location_id, location_name, outcome) |>
  dplyr::mutate(
    RMSE_rank = rank(RMSE, ties.method = "min"),
    MAPE_rank = rank(MAPE, ties.method = "min"),
    mean_rank = (RMSE_rank + MAPE_rank) / 2
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(location_id, outcome, mean_rank, model)

save_table(model_metrics, "Supplement_forecast_model_metrics_all_locations.csv")
save_table(model_ranking, "Supplement_forecast_model_ranking.csv")
save_table(
  dplyr::filter(model_metrics, location_id %in% selected_ids),
  "Table_4_forecast_model_comparison.csv"
)

heatmap_data <- model_metrics |>
  dplyr::filter(
    location_id %in% selected_ids,
    outcome == "Age-standardized incidence rate"
  )
performance_plot <- ggplot2::ggplot(
  heatmap_data, ggplot2::aes(model, location_name, fill = MAPE)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", MAPE)), size = 3) +
  ggplot2::scale_fill_gradient(low = "#EAF2F8", high = "#B03A2E") +
  ggplot2::labs(
    title = "A. Rolling-validation performance", x = NULL, y = NULL,
    fill = "MAPE (%)"
  ) + theme_pcos() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

projection_plot <- future_forecast |>
  dplyr::filter(location_id %in% selected_ids) |>
  ggplot2::ggplot(ggplot2::aes(period_mid, asr, colour = model, group = model)) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::facet_wrap(~location_name, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = FORECAST_PERIOD_MID, labels = FORECAST_PERIOD_LABELS
  ) +
  ggplot2::labs(
    title = "B. Three-model projections, 2024-2038", x = "Projection period",
    y = "Age-standardized incidence rate per 100,000", colour = NULL
  ) + theme_pcos() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 28, hjust = 1),
    legend.position = "bottom"
  )

figure6 <- performance_plot / projection_plot +
  patchwork::plot_layout(heights = c(0.8, 1.7))
save_publication_plot(figure6, "Figure_6_forecast_comparison", width = 13, height = 13)
message("Power5、Poisson和ETS未来预测、滚动验证及Figure 6已完成。")
