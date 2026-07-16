# 年龄标化与ETS训练窗口敏感性分析

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c(
  "data.table", "dplyr", "tidyr", "readr", "ggplot2", "forecast", "patchwork"
))

pcos <- readRDS(PCOS_CLEAN_RDS)
population_historical <- readRDS(POPULATION_RDS)
primary_weights <- readr::read_csv(AGE_WEIGHT_CSV, show_col_types = FALSE) |>
  dplyr::select(age_name, primary_weight = weight)

population_files <- list.files(
  file.path(DATA_DIR, "GBD_population_2023"), pattern = "\\.csv$", full.names = TRUE
)
global_female <- data.table::rbindlist(lapply(population_files, function(path) {
  data.table::fread(
    path,
    select = c("location_name", "sex_name", "age_name", "year", "metric_name", "val"),
    showProgress = FALSE
  )[
    location_name == "Global" & sex_name == "Female" & year == YEAR_END &
      age_name %in% AGE_GROUPS & metric_name == "Number"
  ]
})) |>
  dplyr::distinct(age_name, .keep_all = TRUE) |>
  dplyr::transmute(
    age_name,
    female_population = val,
    female_weight = val / sum(val)
  ) |>
  dplyr::arrange(match(age_name, AGE_GROUPS))
assert_true(nrow(global_female) == length(AGE_GROUPS),
            "女性2023全球人口权重年龄组不完整。")
assert_true(abs(sum(global_female$female_weight) - 1) < 1e-10,
            "女性年龄标准权重之和不等于1。")
save_table(global_female, "Supplement_age_standard_weights_female_2023.csv")

age_rates <- pcos |>
  dplyr::filter(metric_name == "Rate", age_name %in% AGE_GROUPS) |>
  dplyr::select(location_id, location_name, measure, year, age_name, val, lower, upper)

calculate_weighted_rate <- function(data, weights, weight_name, prefix) {
  data |>
    dplyr::inner_join(weights, by = "age_name") |>
    dplyr::group_by(location_id, location_name, measure, year) |>
    dplyr::summarise(
      "{prefix}_asr" := sum(val * .data[[weight_name]]),
      "{prefix}_lower" := sum(lower * .data[[weight_name]]),
      "{prefix}_upper" := sum(upper * .data[[weight_name]]),
      .groups = "drop"
    )
}

primary_asr <- calculate_weighted_rate(
  age_rates, primary_weights, "primary_weight", "primary"
)
female_asr <- calculate_weighted_rate(
  age_rates, global_female, "female_weight", "female"
)
age_standard_sensitivity <- dplyr::left_join(
  primary_asr, female_asr,
  by = c("location_id", "location_name", "measure", "year")
) |>
  dplyr::mutate(
    absolute_difference = female_asr - primary_asr,
    relative_difference_percent =
      100 * (female_asr / pmax(primary_asr, .Machine$double.eps) - 1)
  )

eapc_sensitivity <- age_standard_sensitivity |>
  dplyr::group_by(location_id, location_name, measure) |>
  dplyr::group_modify(function(.x, .y) {
    primary <- eapc_model(dplyr::transmute(.x, year, val = primary_asr))
    female <- eapc_model(dplyr::transmute(.x, year, val = female_asr))
    data.frame(
      primary_EAPC = primary$EAPC,
      primary_EAPC_lower = primary$EAPC_lower,
      primary_EAPC_upper = primary$EAPC_upper,
      female_EAPC = female$EAPC,
      female_EAPC_lower = female$EAPC_lower,
      female_EAPC_upper = female$EAPC_upper,
      EAPC_difference = female$EAPC - primary$EAPC
    )
  }) |>
  dplyr::ungroup()

save_table(
  dplyr::filter(age_standard_sensitivity, year %in% c(YEAR_START, YEAR_END)),
  "Supplement_age_standardization_sensitivity.csv"
)
save_table(eapc_sensitivity, "Supplement_age_standardization_EAPC_sensitivity.csv")

age_standard_summary <- data.frame(
  comparison = "2023 global female weights versus 2023 global both-sex weights",
  maximum_absolute_relative_difference_percent =
    max(abs(age_standard_sensitivity$relative_difference_percent)),
  median_absolute_relative_difference_percent =
    stats::median(abs(age_standard_sensitivity$relative_difference_percent)),
  maximum_absolute_EAPC_difference = max(abs(eapc_sensitivity$EAPC_difference))
)
save_table(age_standard_summary, "QC_07_age_standardization_sensitivity.csv")

incidence_number <- pcos |>
  dplyr::filter(
    measure == "Incidence", metric_name == "Number", age_name %in% AGE_GROUPS
  ) |>
  dplyr::select(location_id, location_name, year, age_name, val)

make_annual_group <- function(ids, group_id, group_name) {
  cases <- incidence_number |>
    dplyr::filter(location_id %in% ids) |>
    dplyr::group_by(year, age_name) |>
    dplyr::summarise(val = sum(val), .groups = "drop")
  population <- population_historical |>
    dplyr::filter(location_id %in% ids) |>
    dplyr::group_by(year, age_name) |>
    dplyr::summarise(population = sum(population), .groups = "drop")
  cases |>
    dplyr::inner_join(population, by = c("year", "age_name")) |>
    dplyr::inner_join(primary_weights, by = "age_name") |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      number = sum(val),
      asr = sum(val / population * RATE_SCALE * primary_weight),
      .groups = "drop"
    ) |>
    dplyr::mutate(location_id = group_id, location_name = group_name)
}

annual_groups <- lapply(seq_len(nrow(COUNTRIES)), function(index) {
  make_annual_group(
    COUNTRIES$location_id[[index]], COUNTRIES$location_id[[index]],
    COUNTRIES$location_name[[index]]
  )
})
annual_groups[[length(annual_groups) + 1L]] <- make_annual_group(
  setdiff(COUNTRIES$location_id, 6L), -1L, "14 land-bordering countries combined"
)
annual_incidence <- dplyr::bind_rows(annual_groups)

fit_ets_window <- function(data, outcome, start_year) {
  training <- data |>
    dplyr::filter(year >= start_year, year <= YEAR_END) |>
    dplyr::arrange(year)
  fit <- forecast::ets(stats::ts(
    training[[outcome]], start = start_year, frequency = 1
  ))
  projected <- forecast::forecast(fit, h = FORECAST_END - YEAR_END, level = 95)
  list(
    annual = data.frame(
      year = (YEAR_END + 1L):FORECAST_END,
      estimate = pmax(as.numeric(projected$mean), 0),
      lower = pmax(as.numeric(projected$lower[, 1]), 0),
      upper = pmax(as.numeric(projected$upper[, 1]), 0)
    ),
    structure = data.frame(
      outcome = outcome,
      training_start = start_year,
      training_end = YEAR_END,
      observations = nrow(training),
      method = fit$method,
      damped = grepl("damped", fit$method, ignore.case = TRUE),
      sigma2 = fit$sigma2,
      AICc = fit$aicc
    )
  )
}

period_starts <- c(2024L, 2029L, 2034L)
period_label_local <- function(year) {
  start <- period_starts[findInterval(year, period_starts)]
  ifelse(year <= 2038L, paste0(start, "-", start + 4L), NA_character_)
}

window_starts <- c(1990L, 2000L, 2010L)
ets_window_results <- list()
ets_structures <- list()
counter <- 1L
for (id in unique(annual_incidence$location_id)) {
  location_data <- dplyr::filter(annual_incidence, location_id == id)
  location_name <- unique(location_data$location_name)[[1]]
  for (start_year in window_starts) {
    number_fit <- fit_ets_window(location_data, "number", start_year)
    asr_fit <- fit_ets_window(location_data, "asr", start_year)
    period_result <- dplyr::left_join(
      number_fit$annual |>
        dplyr::mutate(period = period_label_local(year)) |>
        dplyr::filter(!is.na(period)) |>
        dplyr::group_by(period) |>
        dplyr::summarise(
          number = sum(estimate), number_lower = sum(lower),
          number_upper = sum(upper), .groups = "drop"
        ),
      asr_fit$annual |>
        dplyr::mutate(period = period_label_local(year)) |>
        dplyr::filter(!is.na(period)) |>
        dplyr::group_by(period) |>
        dplyr::summarise(
          asr = mean(estimate), asr_lower = mean(lower),
          asr_upper = mean(upper), .groups = "drop"
        ),
      by = "period"
    ) |>
      dplyr::mutate(
        location_id = id, location_name = location_name,
        training_start = start_year,
        period_mid = as.integer(substr(period, 1, 4)) + 2L
      )
    ets_window_results[[counter]] <- period_result
    ets_structures[[counter]] <- dplyr::bind_rows(
      number_fit$structure, asr_fit$structure
    ) |>
      dplyr::mutate(location_id = id, location_name = location_name)
    counter <- counter + 1L
  }
}

ets_window_results <- dplyr::bind_rows(ets_window_results) |>
  dplyr::arrange(location_id, training_start, period_mid)
ets_structures <- dplyr::bind_rows(ets_structures) |>
  dplyr::arrange(location_id, outcome, training_start)
save_table(ets_window_results, "Supplement_ETS_training_window_sensitivity.csv")
save_table(ets_structures, "Supplement_ETS_model_structures.csv")

recent_trends <- dplyr::bind_rows(lapply(c(2000L, 2010L, 2015L), function(start_year) {
  annual_incidence |>
    dplyr::filter(year >= start_year) |>
    dplyr::group_by(location_id, location_name) |>
    dplyr::group_modify(function(.x, .y) {
      result <- eapc_model(dplyr::transmute(.x, year, val = asr))
      dplyr::mutate(result, trend_start = start_year, trend_end = YEAR_END)
    }) |>
    dplyr::ungroup()
}))
save_table(recent_trends, "Supplement_recent_incidence_EAPC.csv")

primary_projection <- ets_window_results |>
  dplyr::filter(training_start == 1990L, location_id %in% c(6L, -1L))
window_plot <- ets_window_results |>
  dplyr::filter(location_id %in% c(6L, -1L)) |>
  ggplot2::ggplot(ggplot2::aes(
    period_mid, asr, colour = factor(training_start), group = training_start
  )) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::facet_wrap(~location_name, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = c(2026L, 2031L, 2036L),
    labels = c("2024-2028", "2029-2033", "2034-2038")
  ) +
  ggplot2::labs(
    title = "ETS projection sensitivity to the training window",
    subtitle = "Lines show scenarios fitted from 1990, 2000, or 2010 through 2023",
    x = "Projection period", y = "Age-standardized incidence rate per 100,000",
    colour = "Training start"
  ) + theme_pcos() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_publication_plot(
  window_plot, "Figure_S2_ETS_training_window_sensitivity", width = 10, height = 5.8
)

age_plot <- age_standard_sensitivity |>
  dplyr::filter(year == YEAR_END) |>
  ggplot2::ggplot(ggplot2::aes(primary_asr, female_asr, colour = measure)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey45", linetype = 2) +
  ggplot2::geom_point(size = 2, alpha = 0.85) +
  ggplot2::facet_wrap(~measure, scales = "free") +
  ggplot2::labs(
    title = "Sensitivity of 2023 age-standardized rates to standard weights",
    x = "Global both-sex 2023 weights", y = "Global female 2023 weights",
    colour = NULL
  ) + theme_pcos() + ggplot2::theme(legend.position = "none")
save_publication_plot(
  age_plot, "Figure_S1_age_standardization_sensitivity", width = 10, height = 4.8
)

message("年龄标化和ETS训练窗口敏感性分析完成。")
