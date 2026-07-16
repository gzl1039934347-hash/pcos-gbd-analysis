# 1990-2023年PCOS负担人数变化的三因素Shapley分解

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("dplyr", "tidyr", "readr", "ggplot2"))

DECOMPOSITION_REPS <- 1000L

pcos <- readRDS(PCOS_CLEAN_RDS)
population <- readRDS(POPULATION_RDS)

age_numbers <- pcos |>
  dplyr::filter(
    measure %in% MEASURES, metric_name == "Number",
    age_name %in% AGE_GROUPS, year %in% c(YEAR_START, YEAR_END)
  ) |>
  dplyr::select(
    location_id, location_name, measure, year, age_name,
    val, lower, upper
  ) |>
  dplyr::left_join(
    population |>
      dplyr::filter(year %in% c(YEAR_START, YEAR_END)),
    by = c("location_id", "location_name", "year", "age_name")
  )
assert_true(!anyNA(age_numbers), "分解分析输入存在缺失值。")

combined <- age_numbers |>
  dplyr::filter(location_id != 6L) |>
  dplyr::group_by(measure, year, age_name) |>
  dplyr::summarise(
    val = sum(val), lower = sum(lower), upper = sum(upper),
    population = sum(population), .groups = "drop"
  ) |>
  dplyr::mutate(
    location_id = -1L,
    location_name = "14 neighbouring countries combined"
  )
analysis_data <- dplyr::bind_rows(age_numbers, combined)

burden_value <- function(total_population, age_share, age_rate) {
  total_population * sum(age_share * age_rate)
}

shapley_decompose <- function(pop0, pop1, rate0, rate1) {
  n <- c(sum(pop0), sum(pop1))
  share <- list(pop0 / n[[1]], pop1 / n[[2]])
  rate <- list(rate0, rate1)
  factor_names <- c("Population size", "Age composition", "Age-specific rate")
  permutations <- list(
    c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
    c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
  )
  contribution <- setNames(numeric(3), factor_names)
  evaluate <- function(state) {
    burden_value(n[state[[1]] + 1L], share[[state[[2]] + 1L]], rate[[state[[3]] + 1L]])
  }
  for (order in permutations) {
    state <- c(0L, 0L, 0L)
    previous <- evaluate(state)
    for (factor_index in order) {
      state[[factor_index]] <- 1L
      current <- evaluate(state)
      contribution[[factor_index]] <- contribution[[factor_index]] + current - previous
      previous <- current
    }
  }
  contribution <- contribution / length(permutations)
  c(
    contribution,
    `1990 burden` = burden_value(n[[1]], share[[1]], rate[[1]]),
    `2023 burden` = burden_value(n[[2]], share[[2]], rate[[2]])
  )
}

sample_lognormal <- function(val, lower, upper) {
  valid <- val > 0 & lower > 0 & upper > lower
  result <- val
  if (any(valid)) {
    log_se <- (log(upper[valid]) - log(lower[valid])) / (2 * 1.96)
    result[valid] <- stats::rlnorm(
      sum(valid), meanlog = log(val[valid]), sdlog = log_se
    )
  }
  result
}

decomposition_one <- function(data, seed, reps = DECOMPOSITION_REPS) {
  d0 <- data |>
    dplyr::filter(year == YEAR_START) |>
    dplyr::arrange(match(age_name, AGE_GROUPS))
  d1 <- data |>
    dplyr::filter(year == YEAR_END) |>
    dplyr::arrange(match(age_name, AGE_GROUPS))
  point <- shapley_decompose(
    d0$population, d1$population,
    d0$val / d0$population, d1$val / d1$population
  )
  set.seed(seed)
  simulations <- replicate(reps, {
    cases0 <- sample_lognormal(d0$val, d0$lower, d0$upper)
    cases1 <- sample_lognormal(d1$val, d1$lower, d1$upper)
    shapley_decompose(
      d0$population, d1$population,
      cases0 / d0$population, cases1 / d1$population
    )[1:3]
  })
  total_change <- point[["2023 burden"]] - point[["1990 burden"]]
  result <- data.frame(
    component = names(point)[1:3], contribution = point[1:3],
    lower = apply(simulations, 1, stats::quantile, 0.025, na.rm = TRUE),
    upper = apply(simulations, 1, stats::quantile, 0.975, na.rm = TRUE),
    burden_1990 = point[["1990 burden"]],
    burden_2023 = point[["2023 burden"]],
    total_change = total_change,
    stringsAsFactors = FALSE
  )
  result$percentage_of_change <- if (abs(total_change) > 1e-12) {
    100 * result$contribution / total_change
  } else {
    rep(NA_real_, nrow(result))
  }
  result
}

decomposition <- analysis_data |>
  dplyr::group_by(location_id, location_name, measure) |>
  dplyr::group_modify(function(.x, .y) {
    seed <- RANDOM_SEED + abs(as.integer(.y$location_id[[1]])) * 10L +
      match(.y$measure[[1]], MEASURES)
    decomposition_one(.x, seed = seed)
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    component = factor(
      component,
      levels = c("Population size", "Age composition", "Age-specific rate")
    )
  )

reconciliation <- decomposition |>
  dplyr::group_by(location_id, location_name, measure) |>
  dplyr::summarise(
    component_sum = sum(contribution),
    total_change = dplyr::first(total_change),
    relative_error = abs(component_sum - total_change) /
      pmax(abs(total_change), .Machine$double.eps),
    .groups = "drop"
  )
assert_true(max(reconciliation$relative_error) < 1e-8,
            "Shapley分解未能精确还原总负担变化。")
assert_true(all(is.finite(decomposition$contribution)), "分解结果存在非有限值。")
assert_true(all(decomposition$lower <= decomposition$upper), "分解区间顺序错误。")

top3_ids <- age_numbers |>
  dplyr::filter(year == YEAR_END, measure == "Incidence", location_id != 6L) |>
  dplyr::group_by(location_id, location_name) |>
  dplyr::summarise(number = sum(val), .groups = "drop") |>
  dplyr::slice_max(number, n = 3, with_ties = FALSE) |>
  dplyr::pull(location_id)
main_ids <- unique(c(6L, -1L, top3_ids))

save_table(decomposition, "Supplement_decomposition_all_locations.csv")
save_table(reconciliation, "QC_03b_decomposition_reconciliation.csv")
save_table(
  dplyr::filter(decomposition, location_id %in% main_ids),
  "Table_3_decomposition.csv"
)

waterfall_data <- decomposition |>
  dplyr::filter(location_id %in% main_ids) |>
  dplyr::arrange(location_id, measure, component) |>
  dplyr::group_by(location_id, location_name, measure) |>
  dplyr::mutate(
    start = dplyr::lag(cumsum(contribution), default = 0),
    end = cumsum(contribution),
    ymin = pmin(start, end), ymax = pmax(start, end),
    direction = ifelse(contribution >= 0, "Increase", "Decrease")
  ) |>
  dplyr::ungroup()

figure5 <- ggplot2::ggplot(
  waterfall_data,
  ggplot2::aes(x = component, ymin = ymin, ymax = ymax, fill = direction)
) +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = as.numeric(component) - 0.36,
      xmax = as.numeric(component) + 0.36
    ), colour = "white", linewidth = 0.25
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = start + lower, ymax = start + upper),
    width = 0.18, colour = "#303030"
  ) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  ggplot2::facet_wrap(
    ~location_name + measure, ncol = 3, scales = "free_y",
    labeller = ggplot2::labeller(.multi_line = FALSE)
  ) +
  ggplot2::scale_fill_manual(
    values = c(Increase = "#C0392B", Decrease = "#2874A6"), name = NULL
  ) +
  ggplot2::labs(
    title = "Decomposition of changes in PCOS burden, 1990-2023",
    subtitle = "Contributions from population size, age composition, and age-specific rates",
    x = NULL, y = "Change in number"
  ) + theme_pcos() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 28, hjust = 1),
    legend.position = "bottom"
  )
save_publication_plot(figure5, "Figure_5_decomposition", width = 13, height = 13)
message("Shapley负担变化分解及Figure 5已完成。")
