# SII与集中指数分析

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales"))

burden <- readRDS(BURDEN_SUMMARY_RDS)
sdi <- readRDS(SDI_RDS)

analysis_data <- burden |>
  dplyr::filter(year %in% c(YEAR_START, YEAR_END)) |>
  dplyr::select(location_id, location_name, year, measure, val = asr, population) |>
  dplyr::inner_join(dplyr::select(sdi, location_id, year, sdi), by = c("location_id", "year")) |>
  dplyr::group_by(year, measure) |>
  dplyr::mutate(rank = weighted_rank(sdi, population)) |>
  dplyr::ungroup()
assert_true(nrow(analysis_data) == 90L, "不平等分析数据应为90行。")

inequality_results <- analysis_data |>
  dplyr::group_by(year, measure) |>
  dplyr::group_modify(function(.x, .y) {
    fit <- stats::lm(val ~ rank, data = .x, weights = population)
    sii <- stats::coef(fit)[["rank"]]
    sii_ci <- stats::confint(fit, "rank")
    ci <- concentration_index(.x$val, .x$rank, .x$population)
    boot_ci <- bootstrap_inequality(.x)
    data.frame(
      SII = sii,
      SII_lower = sii_ci[1],
      SII_upper = sii_ci[2],
      SII_p = summary(fit)$coefficients["rank", "Pr(>|t|)"],
      concentration_index = ci,
      CI_lower = boot_ci["2.5%", "CI"],
      CI_upper = boot_ci["97.5%", "CI"]
    )
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(SII_p_fdr = stats::p.adjust(SII_p, method = "BH"))
save_table(inequality_results, "Table_2a_SII_concentration_index.csv")

sensitivity <- lapply(c(NA_integer_, KEY_COUNTRIES), function(excluded_id) {
  d <- analysis_data |>
    dplyr::filter(year == YEAR_END)
  label <- "None"
  if (!is.na(excluded_id)) {
    label <- COUNTRIES$location_name[match(excluded_id, COUNTRIES$location_id)]
    d <- dplyr::filter(d, location_id != excluded_id)
  }
  d |>
    dplyr::group_by(measure) |>
    dplyr::group_modify(function(.x, .y) {
      .x$rank <- weighted_rank(.x$sdi, .x$population)
      fit <- stats::lm(val ~ rank, data = .x, weights = population)
      data.frame(
        excluded_country = label,
        SII = stats::coef(fit)[["rank"]],
        concentration_index = concentration_index(.x$val, .x$rank, .x$population)
      )
    }) |>
    dplyr::ungroup()
}) |>
  dplyr::bind_rows()
save_table(sensitivity, "Supplement_inequality_leave_one_out.csv")

curves <- analysis_data |>
  dplyr::group_by(year, measure) |>
  dplyr::group_modify(~ concentration_curve(.x$val, .x$sdi, .x$population)) |>
  dplyr::ungroup()

curve_plot <- ggplot2::ggplot(
  curves,
  ggplot2::aes(population_share, burden_share, colour = factor(year))
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~measure) +
  ggplot2::coord_equal() +
  ggplot2::scale_colour_manual(values = c(`1990` = "#2874A6", `2023` = "#C0392B"), name = "Year") +
  ggplot2::scale_x_continuous(labels = scales::percent_format()) +
  ggplot2::scale_y_continuous(labels = scales::percent_format()) +
  ggplot2::labs(
    title = "A. Concentration curves",
    x = "Cumulative share of women aged 15-49 years",
    y = "Cumulative share of burden"
  ) + theme_pcos()

sii_plot_data <- inequality_results |>
  dplyr::mutate(year = factor(year), measure = factor(measure, levels = MEASURES))
sii_plot <- ggplot2::ggplot(
  sii_plot_data,
  ggplot2::aes(SII, measure, colour = year)
) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = SII_lower, xmax = SII_upper),
    width = 0.16, orientation = "y",
    position = ggplot2::position_dodge(width = 0.45)
  ) +
  ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_colour_manual(values = c(`1990` = "#2874A6", `2023` = "#C0392B"), name = "Year") +
  ggplot2::labs(
    title = "B. Slope index of inequality",
    x = "SII (rate per 100,000)", y = NULL
  ) + theme_pcos()

figure3 <- curve_plot / sii_plot + patchwork::plot_layout(heights = c(1.35, 1))
save_publication_plot(figure3, "Figure_3_inequality", width = 10, height = 10)
message("SII、集中指数及敏感性分析完成。")
