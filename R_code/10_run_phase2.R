# Run Phase 2 scripts in order from the project root.
source("R_code/00_project_helpers.R")
scripts <- c(
  "R_code/03_descriptive_analysis.R",
  "R_code/04_probability_distributions.R",
  "R_code/05_sampling_and_intervals.R",
  "R_code/06_hypothesis_testing.R",
  "R_code/07_correlation_regression.R",
  "R_code/08_regional_analysis.R",
  "R_code/09_time_series_analysis.R"
)
for (s in scripts) {
  cat("\n==== ", s, " ====\n", sep = "")
  sys.source(s, envir = new.env(), chdir = FALSE)
}
cat("\nPhase 2 scripts finished.\n")
