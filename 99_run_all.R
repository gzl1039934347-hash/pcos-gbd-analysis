# 顺序运行PCOS GBD 2023全部分析

scripts <- c(
  "01_clean_pcos.R",
  "02_population_sdi.R",
  "03_descriptive_trends.R",
  "03b_decomposition.R",
  "04_inequality.R",
  "05_frontier.R",
  "06_nordpred.R",
  "07_sensitivity_analysis.R",
  "08_prepare_submission_package.R"
)

for (script in scripts) {
  message("\n========== Running ", script, " ==========")
  source(file.path("E:/GBD_project/scripts", script), encoding = "UTF-8")
}

message("PCOS GBD 2023全流程运行完成。")
