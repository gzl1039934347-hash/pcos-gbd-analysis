# 整理历史人口、SDI、WPP预测人口和2023年年龄权重

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("data.table", "dplyr", "tidyr", "readr", "readxl"))

population_files <- list.files(
  file.path(DATA_DIR, "GBD_population_2023"), pattern = "\\.csv$", full.names = TRUE
)
assert_true(length(population_files) > 0, "未找到GBD 2023人口文件。")

population_all <- data.table::rbindlist(
  lapply(population_files, data.table::fread, showProgress = FALSE),
  use.names = TRUE, fill = TRUE
)

population <- population_all |>
  dplyr::filter(
    location_id %in% COUNTRIES$location_id,
    sex_name == "Female",
    age_name %in% AGE_GROUPS,
    year >= YEAR_START, year <= YEAR_END,
    metric_name == "Number"
  ) |>
  dplyr::transmute(
    location_id,
    location_name = COUNTRIES$location_name[match(location_id, COUNTRIES$location_id)],
    year, age_name, population = val
  ) |>
  dplyr::distinct(location_id, year, age_name, .keep_all = TRUE) |>
  dplyr::arrange(location_id, year, match(age_name, AGE_GROUPS))
assert_true(nrow(population) == 3570L, "历史女性人口应为3,570个组合。")
assert_true(!anyDuplicated(dplyr::select(population, location_id, year, age_name)),
            "历史人口存在重复键。")
saveRDS(population, POPULATION_RDS)
readr::write_csv(population, file.path(DERIVED_DIR, "population_1990_2023.csv"))

global_2023 <- population_all |>
  dplyr::filter(
    location_name == "Global", sex_name == "Both", year == 2023,
    age_name %in% AGE_GROUPS, metric_name == "Number"
  ) |>
  dplyr::transmute(age_name, population = val) |>
  dplyr::distinct(age_name, .keep_all = TRUE)
assert_true(nrow(global_2023) == length(AGE_GROUPS),
            "无法从2023年全球Both-sex人口生成全部15-49岁权重。")
weights <- global_2023 |>
  dplyr::mutate(weight = population / sum(population), weight_percent = 100 * weight) |>
  dplyr::arrange(match(age_name, AGE_GROUPS))
assert_true(abs(sum(weights$weight) - 1) < 1e-10, "年龄权重合计不等于100%。")
readr::write_csv(weights, AGE_WEIGHT_CSV)

sdi_candidates <- list.files(
  file.path(DATA_DIR, "SDI 2023"), pattern = "匹配地区名称后.*\\.xlsx$", full.names = TRUE
)
if (length(sdi_candidates) == 0) {
  sdi_candidates <- list.files(file.path(DATA_DIR, "SDI 2023"), pattern = "^SDI2023\\.xlsx$", full.names = TRUE)
}
assert_true(length(sdi_candidates) >= 1, "未找到SDI 2023工作簿。")
sdi_wide <- readxl::read_excel(sdi_candidates[[1]])
sdi <- sdi_wide |>
  tidyr::pivot_longer(-Location, names_to = "year", values_to = "sdi") |>
  dplyr::mutate(year = as.integer(year)) |>
  dplyr::filter(year >= YEAR_START, year <= YEAR_END) |>
  dplyr::mutate(
    Location = dplyr::recode(
      Location, "Laos" = "Lao People's Democratic Republic",
      "North Korea" = "Democratic People's Republic of Korea"
    )
  ) |>
  dplyr::inner_join(COUNTRIES, by = c("Location" = "location_name")) |>
  dplyr::select(location_id, location_name = Location, iso3, year, sdi) |>
  dplyr::arrange(location_id, year)
assert_true(nrow(sdi) == 510L, "SDI应为15国×34年=510行。")
assert_true(!anyNA(sdi$sdi), "SDI存在缺失值。")
saveRDS(sdi, SDI_RDS)
readr::write_csv(sdi, file.path(DERIVED_DIR, "sdi_1990_2023.csv"))

wpp_path <- file.path(DATA_DIR, "WPP2024_POP_F02_3_POPULATION_5-YEAR_AGE_GROUPS_FEMALE.xlsx")
assert_true(file.exists(wpp_path), "未找到WPP 2024女性5岁年龄组人口文件。")
wpp <- readxl::read_excel(
  wpp_path, sheet = "Medium variant", skip = 16, col_types = "text"
)
wpp_long <- wpp |>
  dplyr::mutate(
    `ISO3 Alpha-code` = as.character(`ISO3 Alpha-code`),
    Year = readr::parse_integer(Year)
  ) |>
  dplyr::filter(
    `ISO3 Alpha-code` %in% COUNTRIES$iso3,
    Year >= YEAR_END + 1L, Year <= FORECAST_END
  ) |>
  dplyr::select(iso3 = `ISO3 Alpha-code`, year = Year, dplyr::all_of(unname(AGE_SHORT))) |>
  tidyr::pivot_longer(dplyr::all_of(unname(AGE_SHORT)), names_to = "age_short", values_to = "population_thousands") |>
  dplyr::mutate(
    population = readr::parse_number(as.character(population_thousands)) * 1000,
    age_name = names(AGE_SHORT)[match(age_short, AGE_SHORT)]
  ) |>
  dplyr::left_join(COUNTRIES, by = "iso3") |>
  dplyr::select(location_id, location_name, iso3, year, age_name, population) |>
  dplyr::arrange(location_id, year, match(age_name, AGE_GROUPS))
assert_true(nrow(wpp_long) == 1785L, "WPP预测人口应为1,785行。")
assert_true(!anyNA(wpp_long$population), "WPP预测人口存在无法解析的值。")
assert_true(!anyDuplicated(dplyr::select(wpp_long, location_id, year, age_name)),
            "WPP预测人口存在重复键。")
saveRDS(wpp_long, WPP_RDS)
readr::write_csv(wpp_long, file.path(DERIVED_DIR, "wpp_female_2024_2040.csv"))

qc <- data.frame(
  check = c("Historical population rows", "SDI rows", "WPP rows", "Age weight sum"),
  value = c(nrow(population), nrow(sdi), nrow(wpp_long), sum(weights$weight_percent)),
  expected = c(3570, 510, 1785, 100)
)
save_table(qc, "QC_02_population_sdi_wpp.csv")
message("人口、SDI、WPP及年龄权重整理完成。")
