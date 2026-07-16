# SDI区域内部低负担基准曲线与发展调整超额YLD差距

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("dplyr", "tidyr", "readr", "ggplot2", "quantreg", "ggrepel"))

burden <- readRDS(BURDEN_SUMMARY_RDS)
sdi <- readRDS(SDI_RDS)
frontier_data <- burden |>
  dplyr::select(location_id, location_name, year, measure, val = asr, population) |>
  dplyr::inner_join(dplyr::select(sdi, location_id, year, sdi), by = c("location_id", "year"))

frontier_2023 <- lapply(MEASURES, function(current_measure) {
  train <- dplyr::filter(frontier_data, measure == current_measure)
  target <- dplyr::filter(train, year == YEAR_END)
  target$frontier <- frontier_predict(train, target$sdi)
  target$efficiency_gap <- pmax(0, target$val - target$frontier)
  target$benchmark_gap <- target$efficiency_gap
  target
}) |>
  dplyr::bind_rows()

set.seed(RANDOM_SEED)
bootstrap_results <- lapply(MEASURES, function(current_measure) {
  train <- dplyr::filter(frontier_data, measure == current_measure)
  target <- dplyr::filter(train, year == YEAR_END)
  country_ids <- unique(train$location_id)
  pred <- replicate(FRONTIER_BOOTSTRAP_REPS, {
    sampled <- sample(country_ids, length(country_ids), replace = TRUE)
    boot <- dplyr::bind_rows(lapply(seq_along(sampled), function(i) {
      x <- train[train$location_id == sampled[[i]], , drop = FALSE]
      x$bootstrap_cluster <- i
      x
    }))
    tryCatch(frontier_predict(boot, target$sdi), error = function(e) rep(NA_real_, nrow(target)))
  })
  target |>
    dplyr::transmute(
      location_id, measure = current_measure,
      frontier_lower = apply(pred, 1, stats::quantile, 0.025, na.rm = TRUE),
      frontier_upper = apply(pred, 1, stats::quantile, 0.975, na.rm = TRUE)
    )
}) |>
  dplyr::bind_rows()

frontier_2023 <- frontier_2023 |>
  dplyr::left_join(bootstrap_results, by = c("location_id", "measure")) |>
  dplyr::mutate(
    gap_lower = pmax(0, val - frontier_upper),
    gap_upper = pmax(0, val - frontier_lower),
    benchmark_gap = efficiency_gap
  ) |>
  dplyr::arrange(measure, dplyr::desc(efficiency_gap))
save_table(frontier_2023, "Table_2b_frontier_efficiency_gap.csv")
save_table(frontier_2023, "Table_2b_regional_benchmark_gap.csv")

leave_one_out <- lapply(KEY_COUNTRIES, function(excluded_id) {
  lapply(MEASURES, function(current_measure) {
    train <- frontier_data |>
      dplyr::filter(measure == current_measure, location_id != excluded_id)
    target <- frontier_data |>
      dplyr::filter(measure == current_measure, year == YEAR_END, location_id != excluded_id)
    target |>
      dplyr::mutate(
        excluded_country = COUNTRIES$location_name[match(excluded_id, COUNTRIES$location_id)],
        frontier = frontier_predict(train, sdi),
        efficiency_gap = pmax(0, val - frontier),
        benchmark_gap = efficiency_gap
      ) |>
      dplyr::select(
        excluded_country, location_id, location_name, measure, sdi, val,
        frontier, efficiency_gap, benchmark_gap
      )
  }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()
save_table(leave_one_out, "Supplement_frontier_leave_one_out.csv")

grid <- data.frame(sdi = seq(min(frontier_data$sdi), max(frontier_data$sdi), length.out = 200))
frontier_lines <- lapply(MEASURES, function(current_measure) {
  train <- dplyr::filter(frontier_data, measure == current_measure)
  data.frame(measure = current_measure, sdi = grid$sdi, frontier = frontier_predict(train, grid$sdi))
}) |>
  dplyr::bind_rows()

figure4 <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = dplyr::filter(frontier_data, year == YEAR_END),
    ggplot2::aes(sdi, val, size = population), alpha = 0.65, colour = "#2874A6"
  ) +
  ggplot2::geom_line(
    data = frontier_lines,
    ggplot2::aes(sdi, frontier), linewidth = 1, colour = "#C0392B"
  ) +
  ggrepel::geom_text_repel(
    data = frontier_2023,
    ggplot2::aes(sdi, val, label = location_name), size = 2.4,
    max.overlaps = 20, box.padding = 0.25
  ) +
  ggplot2::facet_wrap(~measure, scales = "free_y") +
  ggplot2::scale_size_area(
    max_size = 8,
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    name = "Population"
  ) +
  ggplot2::labs(
    title = "Regional development-adjusted lower-burden benchmark, 2023",
    subtitle = "Red line represents the pooled 10th-percentile regional benchmark curve",
    x = "Socio-demographic Index", y = "Age-standardized rate per 100,000"
  ) + theme_pcos()
save_publication_plot(figure4, "Figure_4_frontier", width = 12, height = 6.5)

table2a <- readr::read_csv(file.path(TABLE_DIR, "Table_2a_SII_concentration_index.csv"), show_col_types = FALSE)
table2_combined <- frontier_2023 |>
  dplyr::select(location_id, location_name, measure, sdi, val, frontier,
                frontier_lower, frontier_upper, efficiency_gap, gap_lower, gap_upper) |>
  dplyr::left_join(
    table2a |> dplyr::filter(year == YEAR_END),
    by = "measure"
  )
save_table(table2_combined, "Table_2_inequality_frontier_2023.csv")
message("区域效率前沿、bootstrap及leave-one-country-out分析完成。")
