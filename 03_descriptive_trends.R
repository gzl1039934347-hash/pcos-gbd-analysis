# PCOS负担、年龄模式和EAPC趋势分析

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("dplyr", "tidyr", "readr", "ggplot2", "maps", "patchwork", "scales"))

pcos <- readRDS(PCOS_CLEAN_RDS)
population <- readRDS(POPULATION_RDS)
weights <- readr::read_csv(AGE_WEIGHT_CSV, show_col_types = FALSE)

age_number <- pcos |>
  dplyr::filter(age_name %in% AGE_GROUPS, metric_name == "Number") |>
  dplyr::select(location_id, location_name, year, measure, age_name, val, lower, upper)

age_rate <- pcos |>
  dplyr::filter(age_name %in% AGE_GROUPS, metric_name == "Rate") |>
  dplyr::select(location_id, location_name, year, measure, age_name, val, lower, upper)

numbers <- age_number |>
  dplyr::group_by(location_id, location_name, year, measure) |>
  dplyr::summarise(
    number = sum(val), number_lower = sum(lower), number_upper = sum(upper),
    .groups = "drop"
  )

population_total <- population |>
  dplyr::group_by(location_id, location_name, year) |>
  dplyr::summarise(population = sum(population), .groups = "drop")

asr <- calculate_asr(
  age_rate, weights,
  c("location_id", "location_name", "year", "measure")
) |>
  dplyr::rename(asr = val, asr_lower = lower, asr_upper = upper)

burden <- numbers |>
  dplyr::left_join(population_total, by = c("location_id", "location_name", "year")) |>
  dplyr::mutate(
    crude_rate = number / population * RATE_SCALE,
    crude_lower = number_lower / population * RATE_SCALE,
    crude_upper = number_upper / population * RATE_SCALE
  ) |>
  dplyr::left_join(asr, by = c("location_id", "location_name", "year", "measure")) |>
  dplyr::arrange(measure, location_id, year)
assert_true(nrow(burden) == 1530L, "国家年度负担汇总应为1,530行。")
saveRDS(burden, BURDEN_SUMMARY_RDS)
readr::write_csv(burden, file.path(DERIVED_DIR, "burden_summary.csv"))

eapc <- burden |>
  dplyr::select(location_id, location_name, measure, year, val = asr) |>
  dplyr::group_by(location_id, location_name, measure) |>
  dplyr::group_modify(~ eapc_model(.x)) |>
  dplyr::ungroup()

table1 <- burden |>
  dplyr::filter(year %in% c(YEAR_START, YEAR_END)) |>
  dplyr::select(
    location_id, location_name, measure, year, population,
    number, number_lower, number_upper,
    crude_rate, crude_lower, crude_upper,
    asr, asr_lower, asr_upper
  ) |>
  dplyr::left_join(eapc, by = c("location_id", "location_name", "measure")) |>
  dplyr::arrange(measure, location_id, year)
save_table(table1, "Table_1_PCOS_burden_EAPC.csv")
save_table(eapc, "Supplement_EAPC_1990_2023.csv")

age_pattern <- age_number |>
  dplyr::left_join(population, by = c("location_id", "location_name", "year", "age_name")) |>
  dplyr::mutate(age_band = unname(AGE_BANDS[age_name])) |>
  dplyr::group_by(location_id, location_name, year, measure, age_band) |>
  dplyr::summarise(
    number = sum(val), lower = sum(lower), upper = sum(upper),
    population = sum(population), .groups = "drop"
  ) |>
  dplyr::mutate(
    rate = number / population * RATE_SCALE,
    rate_lower = lower / population * RATE_SCALE,
    rate_upper = upper / population * RATE_SCALE
  )
save_table(age_pattern, "Supplement_age_band_burden.csv")

measure_labels <- c(
  Incidence = "Incidence rate", Prevalence = "Prevalence rate",
  YLD = "YLD rate"
)

figure1 <- ggplot2::ggplot(
  burden,
  ggplot2::aes(year, asr, group = location_name, colour = location_id == 6L)
) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.85) +
  ggplot2::facet_wrap(~measure, scales = "free_y", labeller = ggplot2::as_labeller(measure_labels)) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "#8C9AA6", `TRUE` = "#C0392B"),
    labels = c(`FALSE` = "Neighbouring countries", `TRUE` = "China"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "PCOS burden among women aged 15-49 years, 1990-2023",
    x = NULL, y = "Age-standardized rate per 100,000"
  ) + theme_pcos()
save_publication_plot(figure1, "Figure_1_PCOS_trends", width = 11, height = 6.5)

map_names <- data.frame(
  location_id = c(COUNTRIES$location_id, NA_integer_),
  region = c(
    "China", "Afghanistan", "Bhutan", "India", "Kazakhstan", "Kyrgyzstan",
    "Laos", "Mongolia", "Myanmar", "Nepal", "North Korea", "Pakistan",
    "Russia", "Tajikistan", "Vietnam", "Taiwan"
  ),
  map_label = c(COUNTRIES$location_name, "Taiwan Province, China")
)
map_data <- ggplot2::map_data("world") |>
  dplyr::inner_join(map_names, by = "region") |>
  dplyr::left_join(
    burden |>
      dplyr::filter(year == YEAR_END, measure == "YLD") |>
      dplyr::select(location_id, asr),
    by = "location_id"
  )

map_plot <- ggplot2::ggplot(map_data, ggplot2::aes(long, lat, group = group, fill = asr)) +
  ggplot2::geom_polygon(colour = "white", linewidth = 0.15) +
  ggplot2::geom_polygon(
    data = dplyr::filter(map_data, region == "Taiwan"),
    fill = "#D9D9D9", colour = "#303030", linewidth = 0.35
  ) +
  ggplot2::annotate(
    "segment", x = 121.0, y = 23.6, xend = 128.0, yend = 13.5,
    colour = "#303030", linewidth = 0.3
  ) +
  ggplot2::annotate(
    "text", x = 128.0, y = 12.5,
    label = "Taiwan Province, China\n(no GBD estimate included)",
    hjust = 0.5, vjust = 1, size = 3, family = "Arial", colour = "#303030"
  ) +
  ggplot2::coord_quickmap(xlim = c(20, 150), ylim = c(-5, 60), expand = FALSE) +
  ggplot2::scale_fill_viridis_c(
    option = "C", name = "YLD rate", na.value = "#D9D9D9"
  ) +
  ggplot2::labs(
    title = "A. Age-standardized YLD rate in 2023",
    subtitle = "Taiwan Province, China is shown in grey because no GBD estimate was included"
  ) +
  ggplot2::theme_void(base_size = 11, base_family = "Arial") +
  ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(face = "bold"))

age_plot_data <- age_rate |>
  dplyr::filter(year == YEAR_END, measure == "Incidence") |>
  dplyr::mutate(age_name = factor(age_name, levels = AGE_GROUPS))
age_plot <- ggplot2::ggplot(
  age_plot_data,
  ggplot2::aes(age_name, val, group = location_name, colour = location_id == 6L)
) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.85) +
  ggplot2::geom_point(size = 1) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "#8C9AA6", `TRUE` = "#C0392B"), guide = "none"
  ) +
  ggplot2::labs(
    title = "B. Age-specific incidence rate in 2023",
    x = "Age group", y = "Rate per 100,000"
  ) + theme_pcos() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

figure2 <- map_plot + age_plot + patchwork::plot_layout(widths = c(1.2, 1))
save_publication_plot(figure2, "Figure_2_map_age_pattern", width = 12, height = 5.8)
message("描述性负担、EAPC和Figure 1-2已完成。")
