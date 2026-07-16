# PCOS GBD 2023项目统一配置

PROJECT_ROOT <- "E/"
DATA_DIR <- file.path(PROJECT_ROOT, "data_clean")
SCRIPT_DIR <- file.path(PROJECT_ROOT, "scripts")
OUTPUT_ROOT <- Sys.getenv("PCOS_OUTPUT_ROOT", unset = PROJECT_ROOT)
TABLE_DIR <- file.path(OUTPUT_ROOT, "output_tables")
FIGURE_DIR <- file.path(OUTPUT_ROOT, "output_figures")
DERIVED_DIR <- file.path(DATA_DIR, "pcos_analysis")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)

YEAR_START <- 1990L
YEAR_END <- 2023L
FORECAST_END <- 2040L
RATE_SCALE <- 100000
BOOTSTRAP_REPS <- 1000L
FRONTIER_BOOTSTRAP_REPS <- 300L
RANDOM_SEED <- 20230612L

AGE_GROUPS <- c(
  "15-19 years", "20-24 years", "25-29 years", "30-34 years",
  "35-39 years", "40-44 years", "45-49 years"
)

AGE_SHORT <- setNames(
  c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"),
  AGE_GROUPS
)

AGE_BANDS <- c(
  "15-19 years" = "15-24 years", "20-24 years" = "15-24 years",
  "25-29 years" = "25-34 years", "30-34 years" = "25-34 years",
  "35-39 years" = "35-49 years", "40-44 years" = "35-49 years",
  "45-49 years" = "35-49 years"
)

MEASURES <- c("Incidence", "Prevalence", "YLD")

COUNTRIES <- data.frame(
  location_id = c(6L, 160L, 162L, 163L, 36L, 37L, 12L, 38L, 15L,
                  164L, 7L, 165L, 62L, 39L, 20L),
  location_name = c(
    "China", "Afghanistan", "Bhutan", "India", "Kazakhstan",
    "Kyrgyzstan", "Lao People's Democratic Republic", "Mongolia",
    "Myanmar", "Nepal", "Democratic People's Republic of Korea",
    "Pakistan", "Russian Federation", "Tajikistan", "Viet Nam"
  ),
  iso3 = c("CHN", "AFG", "BTN", "IND", "KAZ", "KGZ", "LAO", "MNG",
           "MMR", "NPL", "PRK", "PAK", "RUS", "TJK", "VNM"),
  stringsAsFactors = FALSE
)

KEY_COUNTRIES <- c(6L, 163L, 62L)

PCOS_CLEAN_RDS <- file.path(DERIVED_DIR, "pcos_clean.rds")
POPULATION_RDS <- file.path(DERIVED_DIR, "population_1990_2023.rds")
SDI_RDS <- file.path(DERIVED_DIR, "sdi_1990_2023.rds")
WPP_RDS <- file.path(DERIVED_DIR, "wpp_female_2024_2040.rds")
AGE_WEIGHT_CSV <- file.path(DERIVED_DIR, "age_weights_15_49_2023.csv")
BURDEN_SUMMARY_RDS <- file.path(DERIVED_DIR, "burden_summary.rds")

source(file.path(SCRIPT_DIR, "utils_pcos.R"), encoding = "UTF-8")
