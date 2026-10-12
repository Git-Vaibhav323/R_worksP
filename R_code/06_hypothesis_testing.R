# Hypothesis tests. H0/H1 stated before each result. alpha = 0.05 unless noted.
# Failure to reject H0 is not proof that H0 is true.
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
subdiv <- dat$subdiv
if (is.character(subdiv$months_complete)) {
  sub_ok <- subdiv[subdiv$months_complete %in% c("TRUE", "true"), ]
} else {
  sub_ok <- subdiv[isTRUE(subdiv$months_complete) | subdiv$months_complete == TRUE, ]
}
sub_ok <- sub_ok[sub_ok$YEAR >= 1901 & sub_ok$YEAR <= 2015, ]

ann <- finite_vec(india$ANNUAL)
jjas <- finite_vec(india$`Jun-Sep`)
ond <- finite_vec(india$`Oct-Dec`)
mam <- finite_vec(india$`Mar-May`)
jf <- finite_vec(india$`Jan-Feb`)
year <- india$YEAR
n <- length(ann)
alpha <- ALPHA

MU0 <- 1180
# Textbook-style all-India annual rainfall of about 118 cm, used as a fixed
# external benchmark, not estimated from this sample.

kerala <- finite_vec(sub_ok$ANNUAL[sub_ok$SUBDIVISION == "Kerala"])
rajas <- finite_vec(sub_ok$ANNUAL[sub_ok$SUBDIVISION == "West Rajasthan"])

t1 <- stats::t.test(ann, mu = MU0, alternative = "two.sided", conf.level = 1 - alpha)
w1 <- stats::wilcox.test(ann, mu = MU0, alternative = "two.sided", exact = FALSE, conf.int = TRUE,
                         conf.level = 1 - alpha)
sw_ann <- stats::shapiro.test(ann)

t2 <- stats::t.test(kerala, rajas, alternative = "two.sided", var.equal = FALSE, conf.level = 1 - alpha)
w2 <- stats::wilcox.test(kerala, rajas, alternative = "two.sided", exact = FALSE, conf.int = TRUE)
sw_k <- stats::shapiro.test(kerala)
sw_r <- stats::shapiro.test(rajas)
var_k <- stats::var(kerala)
var_r <- stats::var(rajas)

t3 <- stats::t.test(jjas, ond, paired = TRUE, alternative = "two.sided", conf.level = 1 - alpha)
w3 <- stats::wilcox.test(jjas, ond, paired = TRUE, alternative = "two.sided", exact = FALSE, conf.int = TRUE)
diff_so <- jjas - ond
sw_diff <- stats::shapiro.test(diff_so)

# Chi-square GOF: which IMD season is the year's wettest? H0: equal probability 1/4.
wettest <- apply(cbind(jf, mam, jjas, ond), 1, which.max)
season_names <- c("Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")
wet_fac <- factor(season_names[wettest], levels = season_names)
wet_tab <- table(wet_fac)
chisq_gof <- stats::chisq.test(wet_tab, p = rep(0.25, 4))

# Chi-square test of independence: early vs late period x monsoon category.
# Categories use the FULL-series JJAS mean so that the cut is not fitted inside each period.
# Years are sequential; treat p-value as approximate.
mid <- 1958
period <- ifelse(year < mid, "1901-1957", "1958-2015")
cat_jjas <- monsoon_category(jjas, mean(jjas))
tab_ind <- table(period, factor(cat_jjas, levels = c("deficient", "normal", "excess")))
chisq_ind <- stats::chisq.test(tab_ind)

# Correlation test (Pearson) between winter and monsoon — not part-whole.
cor_jp <- stats::cor.test(jf, jjas, method = "pearson", conf.level = 1 - alpha)
cor_js <- stats::cor.test(jf, jjas, method = "spearman", exact = FALSE)

# Linear trend coefficient test is reported fully in scripts 07 and 09; H0 stated here.
tr <- stats::cor.test(year, ann, method = "pearson", conf.level = 1 - alpha)

fmt_ci <- function(obj) {
  ci <- obj$conf.int
  if (is.null(ci) || length(ci) < 2) return("NA")
  paste0("[", fmt(ci[1], 3), ", ", fmt(ci[2], 3), "]")
}

tests <- data.frame(
  test_id = 1:8,
  procedure = c(
    "One-sample t-test",
    "Wilcoxon signed-rank (one-sample)",
    "Welch two-sample t-test",
    "Wilcoxon rank-sum",
    "Paired t-test",
    "Wilcoxon signed-rank (paired)",
    "Chi-square goodness-of-fit",
    "Chi-square test of independence"
  ),
  research_question = c(
    "Is mean India annual rainfall different from the 1180 mm textbook benchmark?",
    "Is the median India annual rainfall different from 1180 mm?",
    "Do Kerala and West Rajasthan have different mean annual rainfall (1901-2015 complete years)?",
    "Do Kerala and West Rajasthan differ in annual-rainfall location (stochastic dominance of ranks)?",
    "In the same year, is Jun-Sep rainfall different from Oct-Dec rainfall?",
    "In the same year, is the median of (Jun-Sep minus Oct-Dec) different from zero?",
    "Are the four IMD seasons equally likely to be the wettest season of the year?",
    "Is monsoon deficient/normal/excess independent of early (1901-1957) vs late (1958-2015) period?"
  ),
  H0 = c(
    "mu = 1180 mm",
    "median = 1180 mm",
    "mu_Kerala = mu_WestRajasthan",
    "P(Kerala > Rajasthan) = 0.5 (equal distributions on a shift)",
    "mu_(JJAS - OND) = 0",
    "median(JJAS - OND) = 0",
    "p_JF = p_MAM = p_JJAS = p_OND = 1/4",
    "period and monsoon category are independent"
  ),
  H1 = c(
    "mu != 1180 mm",
    "median != 1180 mm",
    "mu_Kerala != mu_WestRajasthan",
    "the two distributions differ",
    "mu_(JJAS - OND) != 0",
    "median(JJAS - OND) != 0",
    "at least one season probability differs from 1/4",
    "association between period and category"
  ),
  n = c(n, n, paste(length(kerala), length(rajas), sep = "+"),
        paste(length(kerala), length(rajas), sep = "+"),
        n, n, n, n),
  statistic = c(
    unname(t1$statistic), unname(w1$statistic), unname(t2$statistic), unname(w2$statistic),
    unname(t3$statistic), unname(w3$statistic), unname(chisq_gof$statistic), unname(chisq_ind$statistic)
  ),
  df = c(unname(t1$parameter), NA, unname(t2$parameter), NA,
         unname(t3$parameter), NA, unname(chisq_gof$parameter), unname(chisq_ind$parameter)),
  p_value = c(t1$p.value, w1$p.value, t2$p.value, w2$p.value,
              t3$p.value, w3$p.value, chisq_gof$p.value, chisq_ind$p.value),
  ci_95 = c(fmt_ci(t1), fmt_ci(w1), fmt_ci(t2), "see output file",
            fmt_ci(t3), fmt_ci(w3), NA, NA),
  alpha = alpha,
  decision = c(
    ifelse(t1$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(w1$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(t2$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(w2$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(t3$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(w3$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(chisq_gof$p.value < alpha, "Reject H0", "Do not reject H0"),
    ifelse(chisq_ind$p.value < alpha, "Reject H0", "Do not reject H0")
  ),
  stringsAsFactors = FALSE
)
write.csv(tests, file.path("results", "tables", "hypothesis_test_results.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(tab_ind),
          file.path("results", "tables", "chisq_period_by_monsoon_category.csv"))
write.csv(as.data.frame(wet_tab), file.path("results", "tables", "wettest_season_counts.csv"),
          row.names = FALSE)

cor_tbl <- data.frame(
  pair = c("Jan-Feb vs Jun-Sep (Pearson)", "Jan-Feb vs Jun-Sep (Spearman)",
           "YEAR vs annual (Pearson trend association)"),
  estimate = c(unname(cor_jp$estimate), unname(cor_js$estimate), unname(tr$estimate)),
  p_value = c(cor_jp$p.value, cor_js$p.value, tr$p.value),
  ci_lower = c(cor_jp$conf.int[1], NA, tr$conf.int[1]),
  ci_upper = c(cor_jp$conf.int[2], NA, tr$conf.int[2]),
  stringsAsFactors = FALSE
)
write.csv(cor_tbl, file.path("results", "tables", "correlation_significance_tests.csv"), row.names = FALSE)

lines <- c(
  "HYPOTHESIS TESTING REPORT",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Significance level alpha =", alpha, "for all tests in this file."),
  "Not rejecting H0 does not prove H0. P-values are not posterior probabilities.",
  "",
  "================================================================",
  "TEST 1. One-sample t-test",
  "================================================================",
  "Research question: Is the mean of India area-weighted annual rainfall different",
  "from 1180 mm (a commonly quoted ~118 cm all-India annual figure used as a fixed benchmark)?",
  "H0: mu = 1180 mm.  H1: mu != 1180 mm.",
  paste("n =", n, "years (1901-2015). Sample mean =", fmt(mean(ann), 2), "mm."),
  "Assumptions: observations iid normal. Shapiro-Wilk on annual totals:",
  paste0("  W = ", fmt(as.numeric(sw_ann$statistic), 4), ", p = ", fmt_p(sw_ann$p.value), "."),
  "Mild skewness is tolerated because n = 115 and the CLT supports the t approximation.",
  "Years are a time series, so iid is an approximation (see script 09).",
  paste0("t = ", fmt(as.numeric(t1$statistic), 3), ", df = ", fmt(as.numeric(t1$parameter), 1),
         ", p = ", fmt_p(t1$p.value), "."),
  paste("95% CI for mu:", fmt_ci(t1), "mm."),
  paste("Decision:", ifelse(t1$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  "Interpretation: compare the interval to 1180 mm; a rejection would mean the historical",
  "mean of THIS series is not 1180 mm, not that textbooks are 'wrong' in general.",
  "",
  "TEST 1b. Wilcoxon signed-rank (robust location)",
  paste0("V = ", fmt(as.numeric(w1$statistic), 1), ", p = ", fmt_p(w1$p.value),
         ", 95% CI (pseudomedian) ", fmt_ci(w1), "."),
  paste("Decision:", ifelse(w1$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  "",
  "================================================================",
  "TEST 2. Welch two-sample t-test (independent regions)",
  "================================================================",
  "Research question: Do Kerala and West Rajasthan differ in mean annual rainfall?",
  "H0: mu_K = mu_R.  H1: mu_K != mu_R.",
  paste("n_Kerala =", length(kerala), "; n_West Rajasthan =", length(rajas),
        "; complete years 1901-2015."),
  "These are different places, so a two-sample design is appropriate. Do NOT use this",
  "design on two seasons from the same year (that would be paired; see Test 3).",
  paste0("Shapiro-Wilk Kerala p = ", fmt_p(sw_k$p.value),
         "; West Rajasthan p = ", fmt_p(sw_r$p.value), "."),
  paste0("Sample variances: Kerala ", fmt(var_k, 1), ", West Rajasthan ", fmt(var_r, 1),
         " (unequal; Welch is used, not Student's equal-variance t)."),
  paste0("t = ", fmt(as.numeric(t2$statistic), 3), ", df = ", fmt(as.numeric(t2$parameter), 2),
         ", p = ", fmt_p(t2$p.value), "."),
  paste("95% CI for mu_K - mu_R:", fmt_ci(t2), "mm."),
  paste("Decision:", ifelse(t2$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  "Interpretation: a large difference is expected climatologically (wet coast vs arid NW).",
  "The test documents that difference; it is not a surprising scientific discovery.",
  paste0("Wilcoxon rank-sum W = ", fmt(as.numeric(w2$statistic), 1), ", p = ", fmt_p(w2$p.value), "."),
  "",
  "================================================================",
  "TEST 3. Paired t-test (same year, two seasons)",
  "================================================================",
  "Research question: Within a year, does southwest-monsoon rainfall differ from post-monsoon?",
  "H0: mean(JJAS - OND) = 0.  H1: not 0.",
  paste("n =", n, "paired years. Pairing is required because both seasons belong to the same year."),
  paste0("Shapiro-Wilk on differences: W = ", fmt(as.numeric(sw_diff$statistic), 4),
         ", p = ", fmt_p(sw_diff$p.value), "."),
  paste0("t = ", fmt(as.numeric(t3$statistic), 3), ", df = ", fmt(as.numeric(t3$parameter), 1),
         ", p = ", fmt_p(t3$p.value), "."),
  paste("95% CI for mean seasonal difference:", fmt_ci(t3), "mm."),
  paste("Decision:", ifelse(t3$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  paste0("Paired Wilcoxon p = ", fmt_p(w3$p.value), "."),
  "",
  "================================================================",
  "TEST 4. Chi-square goodness-of-fit (categorical frequencies, not raw mm)",
  "================================================================",
  "Each year is classified by which of the four IMD seasons has the largest total.",
  "H0: each season is equally likely to be wettest (probability 1/4).",
  "This is a valid categorical GOF. Chi-square is NOT applied to raw millimetres.",
  paste("Counts:", paste(paste(names(wet_tab), as.integer(wet_tab), sep = "="), collapse = ", "), "."),
  paste0("X^2 = ", fmt(as.numeric(chisq_gof$statistic), 2), ", df = ", chisq_gof$parameter,
         ", p = ", fmt_p(chisq_gof$p.value), "."),
  paste("Expected count under H0 =", n / 4, "per season. Decision:",
        ifelse(chisq_gof$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  "A rejection is climatologically expected because the southwest monsoon dominates India.",
  "",
  "================================================================",
  "TEST 5. Chi-square independence (period x monsoon category)",
  "================================================================",
  "H0: early/late period is independent of deficient/normal/excess JJAS class.",
  "Classes use 90% and 110% of the 1901-2015 JJAS mean.",
  "Expected cell counts:",
  paste(capture.output(print(round(chisq_ind$expected, 2))), collapse = "\n"),
  paste0("X^2 = ", fmt(as.numeric(chisq_ind$statistic), 2), ", df = ", chisq_ind$parameter,
         ", p = ", fmt_p(chisq_ind$p.value), "."),
  paste("Decision:", ifelse(chisq_ind$p.value < alpha, "Reject H0.", "Do not reject H0.")),
  "Limitation: adjacent years are dependent, so the chi-square nominal p-value is approximate.",
  "This is not a test of climate change causation.",
  "",
  "Correlation tests (detail in script 07):",
  paste0("Pearson JF vs JJAS r = ", fmt(unname(cor_jp$estimate), 3), ", p = ", fmt_p(cor_jp$p.value), "."),
  paste0("Pearson YEAR vs annual r = ", fmt(unname(tr$estimate), 3), ", p = ", fmt_p(tr$p.value), "."),
  "",
  "One-way ANOVA and post-hoc regional comparisons are in 08_regional_analysis.R",
  "because they require the subdivision panel and assumption checks by region.",
  "Regression coefficient tests are in 07_correlation_regression.R."
)
write_txt(file.path("results", "statistical_outputs", "hypothesis_testing_report.txt"), lines)

cat("06_hypothesis_testing.R complete.\n")
