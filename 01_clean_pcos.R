# 合并并清洗PCOS GBD 2023数据

source("E:/GBD_project/scripts/00_config.R", encoding = "UTF-8")
required_packages(c("data.table", "dplyr", "readr"))

all_csv <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)

is_pcos_file <- function(path) {
  sample_row <- tryCatch(
    data.table::fread(path, nrows = 1, showProgress = FALSE),
    error = function(e) NULL
  )
  !is.null(sample_row) && "cause_name" %in% names(sample_row) &&
    identical(sample_row$cause_name[[1]], "Polycystic ovarian syndrome")
}

pcos_files <- all_csv[vapply(all_csv, is_pcos_file, logical(1))]
assert_true(length(pcos_files) > 0, "未找到PCOS CSV文件。")
message("读取PCOS文件数: ", length(pcos_files))

raw <- data.table::rbindlist(
  lapply(pcos_files, data.table::fread, showProgress = FALSE),
  use.names = TRUE, fill = TRUE, idcol = "source_file"
)
raw$source_file <- basename(pcos_files[raw$source_file])

country_ids <- COUNTRIES$location_id
source_check <- raw |>
  dplyr::filter(
    cause_name == "Polycystic ovarian syndrome",
    location_id %in% country_ids,
    sex_name == "Female",
    year >= YEAR_START, year <= YEAR_END,
    age_name %in% AGE_GROUPS,
    metric_name %in% c("Number", "Rate")
  ) |>
  dplyr::mutate(measure = normalize_measure(measure_name)) |>
  dplyr::filter(measure %in% c(MEASURES, "DALY"))

source_key <- c(
  "location_id", "sex_id", "age_id", "cause_id", "measure_id",
  "metric_id", "year"
)
source_unique <- source_check |> dplyr::distinct(dplyr::across(dplyr::all_of(source_key)))
assert_true(
  nrow(source_unique) == 28560L,
  paste0("来源数据组合应为28,560，实际为", nrow(source_unique), "。")
)

clean <- raw |>
  dplyr::filter(
    cause_name == "Polycystic ovarian syndrome",
    location_id %in% country_ids,
    sex_name == "Female",
    year >= YEAR_START, year <= YEAR_END,
    age_name %in% c(AGE_GROUPS, "15-49 years", "Age-standardized"),
    metric_name %in% c("Number", "Rate")
  ) |>
  dplyr::mutate(measure = normalize_measure(measure_name)) |>
  dplyr::filter(measure %in% MEASURES) |>
  dplyr::mutate(
    location_name = COUNTRIES$location_name[match(location_id, COUNTRIES$location_id)]
  ) |>
  dplyr::select(
    source_file, measure_id, measure, location_id, location_name,
    sex_id, sex_name, age_id, age_name, cause_id, cause_name,
    metric_id, metric_name, year, val, lower, upper
  )

analysis_key <- c(
  "measure", "location_id", "sex_id", "age_id", "cause_id",
  "metric_id", "year"
)

conflicts <- clean |>
  dplyr::group_by(dplyr::across(dplyr::all_of(analysis_key))) |>
  dplyr::summarise(
    n_values = dplyr::n_distinct(paste(val, lower, upper)),
    .groups = "drop"
  ) |>
  dplyr::filter(n_values > 1)
assert_true(nrow(conflicts) == 0, "发现同一维度键对应不同估计值，不能自动去重。")

duplicate_rows <- nrow(clean) - nrow(dplyr::distinct(clean, dplyr::across(dplyr::all_of(analysis_key))))
clean <- clean |>
  dplyr::arrange(location_id, measure, year, age_id, metric_id) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(analysis_key)), .keep_all = TRUE)

final_age_specific <- clean |>
  dplyr::filter(age_name %in% AGE_GROUPS) |>
  dplyr::distinct(dplyr::across(dplyr::all_of(analysis_key)))
assert_true(
  nrow(final_age_specific) == 21420L,
  paste0("最终PCOS年龄别组合应为21,420，实际为", nrow(final_age_specific), "。")
)
assert_true(all(!is.na(clean$val) & !is.na(clean$lower) & !is.na(clean$upper)),
            "PCOS数据存在缺失估计值或不确定性区间。")

saveRDS(clean, PCOS_CLEAN_RDS)
readr::write_csv(clean, file.path(DERIVED_DIR, "pcos_clean.csv"))

qc <- data.frame(
  check = c(
    "PCOS source files", "Source four-measure age-specific keys",
    "Retained three-measure age-specific keys", "Duplicate rows removed",
    "Conflicting duplicate keys", "Missing val/lower/upper"
  ),
  value = c(length(pcos_files), nrow(source_unique), nrow(final_age_specific),
            duplicate_rows, nrow(conflicts),
            sum(is.na(dplyr::select(clean, val, lower, upper)))),
  expected = c(NA, 28560, 21420, NA, 0, 0)
)
save_table(qc, "QC_01_pcos_cleaning.csv")
message("PCOS清洗完成，去除重复行: ", duplicate_rows)
