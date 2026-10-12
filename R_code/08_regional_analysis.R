# Regional summaries and group comparisons. Common period 1901-2015, complete years.
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
subdiv <- dat$subdiv
if (is.character(subdiv$months_complete)) {
  sub_ok <- subdiv[subdiv$months_complete %in% c("TRUE", "true"), ]
} else {
  sub_ok <- subdiv[subdiv$months_complete == TRUE, ]
}
sub_ok <- sub_ok[sub_ok$YEAR >= 1901 & sub_ok$YEAR <= 2015, ]

regs <- sort(unique(sub_ok$SUBDIVISION))
reg_rows <- lapply(regs, function(rg) {
  d <- sub_ok[sub_ok$SUBDIVISION == rg, ]
  a <- desc_stats(d$ANNUAL, rg)
  j <- desc_stats(d$JJAS, rg)
  data.frame(
    SUBDIVISION = rg,
    n_years = a$n,
    year_min = min(d$YEAR),
    year_max = max(d$YEAR),
    annual_mean = a$mean,
    annual_sd = a$sd,
    annual_cv_pct = a$cv_pct,
    annual_median = a$median,
    annual_min = a$min,
    annual_max = a$max,
    annual_skewness = a$skewness,
    jjas_mean = j$mean,
    jjas_sd = j$sd,
    jjas_cv_pct = j$cv_pct,
    stringsAsFactors = FALSE
  )
})
reg_tbl <- do.call(rbind, reg_rows)
reg_tbl <- reg_tbl[order(-reg_tbl$annual_mean), ]
write.csv(reg_tbl, file.path("results", "tables", "regional_annual_summaries.csv"), row.names = FALSE)

# Monthly regional means (complete years)
month_means <- aggregate(sub_ok[MONTH_COLS],
                         by = list(SUBDIVISION = sub_ok$SUBDIVISION),
                         FUN = mean, na.rm = TRUE)
write.csv(month_means, file.path("results", "tables", "regional_monthly_means.csv"), row.names = FALSE)

# Selected regions spanning wet coast, interior, arid NW, and SE (NE monsoon).
selected <- c("Kerala", "Coastal Karnataka", "Gangetic West Bengal",
              "Tamil Nadu", "Punjab", "West Rajasthan")
sel <- sub_ok[sub_ok$SUBDIVISION %in% selected, ]
sel$SUBDIVISION <- factor(sel$SUBDIVISION, levels = selected)

# Balanced-ish: drop years not present in all six (inner join on YEAR)
year_counts <- table(sel$YEAR)
common_years <- as.integer(names(year_counts)[year_counts == length(selected)])
sel_bal <- sel[sel$YEAR %in% common_years, ]

aov_fit <- stats::aov(ANNUAL ~ SUBDIVISION, data = sel_bal)
aov_s <- summary(aov_fit)[[1]]
welch <- stats::oneway.test(ANNUAL ~ SUBDIVISION, data = sel_bal, var.equal = FALSE)
kw <- stats::kruskal.test(ANNUAL ~ SUBDIVISION, data = sel_bal)
flig <- stats::fligner.test(ANNUAL ~ SUBDIVISION, data = sel_bal)

# Shapiro by group
shap <- sapply(selected, function(rg) {
  x <- sel_bal$ANNUAL[sel_bal$SUBDIVISION == rg]
  stats::shapiro.test(x)$p.value
})

# Tukey HSD (assumes equal variance; reported with a caveat if Fligner rejects)
tuk <- stats::TukeyHSD(aov_fit, conf.level = 0.95)
tuk_df <- as.data.frame(tuk$SUBDIVISION)
tuk_df$pair <- rownames(tuk_df)
write.csv(tuk_df, file.path("results", "tables", "anova_tukey_hsd.csv"), row.names = FALSE)

# Pairwise Welch t-tests with Bonferroni (more appropriate if variances differ)
pw <- stats::pairwise.t.test(sel_bal$ANNUAL, sel_bal$SUBDIVISION,
                             pool.sd = FALSE, p.adjust.method = "bonferroni")
pw_df <- as.data.frame(as.table(pw$p.value))
names(pw_df) <- c("group1", "group2", "p_bonferroni")
pw_df <- pw_df[!is.na(pw_df$p_bonferroni), ]
write.csv(pw_df, file.path("results", "tables", "anova_pairwise_welch_bonferroni.csv"), row.names = FALSE)

# Pairwise Wilcoxon with Bonferroni
pw_w <- stats::pairwise.wilcox.test(sel_bal$ANNUAL, sel_bal$SUBDIVISION,
                                    p.adjust.method = "bonferroni", exact = FALSE)
pw_w_df <- as.data.frame(as.table(pw_w$p.value))
names(pw_w_df) <- c("group1", "group2", "p_bonferroni")
pw_w_df <- pw_w_df[!is.na(pw_w_df$p_bonferroni), ]
write.csv(pw_w_df, file.path("results", "tables", "anova_pairwise_wilcoxon_bonferroni.csv"),
          row.names = FALSE)

anova_sum <- data.frame(
  test = c("Classical one-way ANOVA", "Welch one-way ANOVA", "Kruskal-Wallis", "Fligner-Killeen (variances)"),
  statistic = c(aov_s[1, "F value"], unname(welch$statistic), unname(kw$statistic), unname(flig$statistic)),
  df1 = c(aov_s[1, "Df"], welch$parameter[1], NA, flig$parameter),
  df2 = c(aov_s[2, "Df"], welch$parameter[2], kw$parameter, NA),
  p_value = c(aov_s[1, "Pr(>F)"], welch$p.value, kw$p.value, flig$p.value),
  n_years_balanced = length(common_years),
  n_regions = length(selected),
  stringsAsFactors = FALSE
)
write.csv(anova_sum, file.path("results", "tables", "anova_regional_annual.csv"), row.names = FALSE)

# Regional linear slopes (descriptive) on the balanced window
slope_rows <- lapply(selected, function(rg) {
  d <- sel_bal[sel_bal$SUBDIVISION == rg, ]
  m <- stats::lm(ANNUAL ~ YEAR, data = d)
  s <- summary(m)
  data.frame(
    SUBDIVISION = rg,
    n = nrow(d),
    slope_mm_per_year = unname(coef(m)[2]),
    slope_se = s$coefficients[2, 2],
    slope_p = s$coefficients[2, 4],
    r_squared = s$r.squared,
    stringsAsFactors = FALSE
  )
})
write.csv(do.call(rbind, slope_rows),
          file.path("results", "tables", "regional_ols_trends.csv"), row.names = FALSE)

write_txt(file.path("results", "statistical_outputs", "regional_anova_notes.txt"), c(
  "REGIONAL ANALYSIS AND ANOVA",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "SCOPE",
  "All regional means use complete-month years in 1901-2015 so that annual totals are defined",
  "and the window matches the national series. Subdivisions have different rainfall climates;",
  "a larger mean in Kerala than in West Rajasthan is a climatological fact, not a data error.",
  "Equal-weighted averages of subdivisions are NOT the national area-weighted series.",
  "",
  paste("Selected regions:", paste(selected, collapse = ", ")),
  paste("Balanced common years with all six present:", length(common_years),
        "(", min(common_years), "-", max(common_years), ")."),
  "",
  "ASSUMPTIONS",
  "Classical ANOVA needs independent groups, similar variances, and approximate normality.",
  "Regions are different places (independent groups is reasonable). Years within a region",
  "are a time series (independence across years is approximate).",
  paste0("Fligner-Killeen variance test: chi-square = ", fmt(unname(flig$statistic), 2),
         ", p = ", fmt_p(flig$p.value), "."),
  paste("Shapiro-Wilk p-values by region:",
        paste(paste(names(shap), fmt_p(shap), sep = "="), collapse = "; "), "."),
  "If Fligner rejects equal variances, prefer Welch ANOVA and Bonferroni pairwise Welch tests.",
  "Kruskal-Wallis does not assume normality or equal variances (still assumes similar shapes).",
  "",
  "RESULTS",
  paste0("Classical ANOVA F = ", fmt(aov_s[1, "F value"], 2),
         ", df = ", aov_s[1, "Df"], ", ", aov_s[2, "Df"],
         ", p = ", fmt_p(aov_s[1, "Pr(>F)"]), "."),
  paste0("Welch ANOVA F = ", fmt(unname(welch$statistic), 2),
         ", p = ", fmt_p(welch$p.value), "."),
  paste0("Kruskal-Wallis chi-square = ", fmt(unname(kw$statistic), 2),
         ", df = ", kw$parameter, ", p = ", fmt_p(kw$p.value), "."),
  paste("Decision at alpha = 0.05:",
        ifelse(welch$p.value < ALPHA, "Reject equal regional means.", "Do not reject equal means.")),
  "Post-hoc: Tukey HSD (equal-variance) and Bonferroni-adjusted pairwise Welch and Wilcoxon",
  "are saved as CSV. Interpret Tukey cautiously if variances differ.",
  "",
  "Regional OLS slopes vs year are descriptive; Mann-Kendall in script 09 is preferred",
  "when residual autocorrelation is a concern. No causal climate attribution is made."
))

ord <- order(reg_tbl$annual_mean)
save_png("regional_mean_annual_rainfall.png", {
  graphics::par(mar = c(10, 5, 3, 1))
  graphics::barplot(reg_tbl$annual_mean[ord], names.arg = reg_tbl$SUBDIVISION[ord],
                    las = 2, col = "#2B8CBE", border = NA, cex.names = 0.55,
                    ylab = "Mean annual rainfall (mm)",
                    main = "Mean annual rainfall by IMD subdivision (complete years, 1901-2015)")
})

save_png("regional_cv_annual_rainfall.png", {
  graphics::par(mar = c(10, 5, 3, 1))
  cvord <- order(reg_tbl$annual_cv_pct)
  graphics::barplot(reg_tbl$annual_cv_pct[cvord], names.arg = reg_tbl$SUBDIVISION[cvord],
                    las = 2, col = "#E34A33", border = NA, cex.names = 0.55,
                    ylab = "Coefficient of variation (%)",
                    main = "Relative variability of annual rainfall by subdivision, 1901-2015")
})

save_png("selected_regions_annual_boxplots.png", {
  graphics::boxplot(ANNUAL ~ SUBDIVISION, data = sel_bal, las = 2, col = "#A6BDDB",
                    ylab = "Annual rainfall (mm)", xlab = "",
                    main = "Annual rainfall in six selected subdivisions (common years)")
})

cat("08_regional_analysis.R complete.\n")
