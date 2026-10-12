# Correlation (not causation) and regression. No part-whole pairs (annual vs months).
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
seas <- india[, c("Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")]
names(seas) <- c("JanFeb", "MarMay", "JunSep", "OctDec")
ann <- india$ANNUAL
year <- india$YEAR

pear <- stats::cor(seas, method = "pearson")
spear <- stats::cor(seas, method = "spearman")
cvm <- stats::cov(seas)

write.csv(as.data.frame(pear), file.path("results", "tables", "seasonal_pearson_correlation.csv"))
write.csv(as.data.frame(spear), file.path("results", "tables", "seasonal_spearman_correlation.csv"))
write.csv(as.data.frame(cvm), file.path("results", "tables", "seasonal_covariance.csv"))

pairs <- list(
  c("JanFeb", "MarMay"), c("JanFeb", "JunSep"), c("JanFeb", "OctDec"),
  c("MarMay", "JunSep"), c("MarMay", "OctDec"), c("JunSep", "OctDec")
)
pair_rows <- lapply(pairs, function(pr) {
  x <- seas[[pr[1]]]; y <- seas[[pr[2]]]
  pe <- stats::cor.test(x, y, method = "pearson")
  sp <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
  data.frame(
    var_x = pr[1], var_y = pr[2],
    pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
    pearson_ci_lower = pe$conf.int[1], pearson_ci_upper = pe$conf.int[2],
    spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
    n = length(x),
    stringsAsFactors = FALSE
  )
})
pair_tbl <- do.call(rbind, pair_rows)
write.csv(pair_tbl, file.path("results", "tables", "seasonal_pairwise_correlations.csv"), row.names = FALSE)

# Simple linear regression: descriptive trend of annual rainfall on year.
m1 <- stats::lm(ANNUAL ~ YEAR, data = india)
s1 <- summary(m1)
ci1 <- stats::confint(m1, level = 0.95)
e1 <- stats::residuals(m1)
f1 <- stats::fitted(m1)
dw1 <- durbin_watson(e1)
lag1 <- stats::cor.test(e1[-length(e1)], e1[-1])
cook1 <- stats::cooks.distance(m1)
infl_n <- sum(cook1 > 4 / length(cook1), na.rm = TRUE)
bp <- stats::lm(e1^2 ~ f1)
bp_p <- summary(bp)$coefficients[2, 4]
sw_e <- stats::shapiro.test(e1)

# Multiple regression with a scientific rationale: monsoon vs year and winter.
# Winter is a different season, not a component of Jun-Sep.
m2 <- stats::lm(`Jun-Sep` ~ YEAR + `Jan-Feb`, data = india)
s2 <- summary(m2)
ci2 <- stats::confint(m2, level = 0.95)
e2 <- stats::residuals(m2)
dw2 <- durbin_watson(e2)
vif_approx <- {
  r <- summary(stats::lm(YEAR ~ `Jan-Feb`, data = india))$r.squared
  c(YEAR = 1 / (1 - r), JanFeb = 1 / (1 - r))
}

coef_tbl <- rbind(
  data.frame(model = "annual ~ year",
             term = rownames(s1$coefficients),
             estimate = s1$coefficients[, 1],
             std_error = s1$coefficients[, 2],
             t_value = s1$coefficients[, 3],
             p_value = s1$coefficients[, 4],
             ci95_lower = ci1[, 1],
             ci95_upper = ci1[, 2],
             stringsAsFactors = FALSE),
  data.frame(model = "JunSep ~ year + JanFeb",
             term = rownames(s2$coefficients),
             estimate = s2$coefficients[, 1],
             std_error = s2$coefficients[, 2],
             t_value = s2$coefficients[, 3],
             p_value = s2$coefficients[, 4],
             ci95_lower = ci2[, 1],
             ci95_upper = ci2[, 2],
             stringsAsFactors = FALSE)
)
write.csv(coef_tbl, file.path("results", "tables", "regression_coefficients.csv"), row.names = FALSE)

fit_tbl <- data.frame(
  model = c("annual ~ year", "JunSep ~ year + JanFeb"),
  n = c(nobs(m1), nobs(m2)),
  r_squared = c(s1$r.squared, s2$r.squared),
  adj_r_squared = c(s1$adj.r.squared, s2$adj.r.squared),
  residual_se = c(s1$sigma, s2$sigma),
  F_stat = c(s1$fstatistic[1], s2$fstatistic[1]),
  F_df1 = c(s1$fstatistic[2], s2$fstatistic[2]),
  F_df2 = c(s1$fstatistic[3], s2$fstatistic[3]),
  F_p = c(
    stats::pf(s1$fstatistic[1], s1$fstatistic[2], s1$fstatistic[3], lower.tail = FALSE),
    stats::pf(s2$fstatistic[1], s2$fstatistic[2], s2$fstatistic[3], lower.tail = FALSE)
  ),
  durbin_watson = c(dw1, dw2),
  residual_lag1_p = c(lag1$p.value, NA_real_),
  shapiro_resid_p = c(sw_e$p.value, stats::shapiro.test(e2)$p.value),
  hetero_p_squared_resid_on_fitted = c(bp_p, summary(stats::lm(e2^2 ~ stats::fitted(m2)))$coefficients[2, 4]),
  n_cook_gt_4_over_n = c(infl_n, sum(stats::cooks.distance(m2) > 4 / nobs(m2))),
  stringsAsFactors = FALSE
)
write.csv(fit_tbl, file.path("results", "tables", "regression_model_fit.csv"), row.names = FALSE)

cook_tbl <- data.frame(
  YEAR = india$YEAR,
  ANNUAL = india$ANNUAL,
  residual = e1,
  fitted = f1,
  cook_d = as.numeric(cook1),
  influential_4_over_n = as.numeric(cook1) > 4 / length(cook1)
)
write.csv(cook_tbl, file.path("results", "tables", "regression_influence_annual_on_year.csv"),
          row.names = FALSE)

write_txt(file.path("results", "statistical_outputs", "correlation_regression_notes.txt"), c(
  "CORRELATION AND REGRESSION",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "CORRELATION RULES USED",
  "Pearson measures linear association; Spearman measures monotonic association.",
  "Neither implies causation. Seasonal totals JF, MAM, JJAS, OND are different seasons,",
  "so their correlations are not part-whole artefacts. Annual vs Jun-Sep is omitted",
  "because Jun-Sep is a large numerical component of annual rainfall.",
  "",
  paste0("Largest |Pearson| among seasons is documented in seasonal_pairwise_correlations.csv."),
  paste0("Example: Jan-Feb vs Jun-Sep Pearson r = ",
         fmt(pair_tbl$pearson_r[pair_tbl$var_x == "JanFeb" & pair_tbl$var_y == "JunSep"], 3),
         " (p = ",
         fmt_p(pair_tbl$pearson_p[pair_tbl$var_x == "JanFeb" & pair_tbl$var_y == "JunSep"]),
         ")."),
  "A weak winter-monsoon correlation, if present, is association only.",
  "Outliers and nonlinearity: see scatterplots. Time dependence can inflate significance.",
  "",
  "MODEL 1  ANNUAL = beta0 + beta1 * YEAR + error",
  paste0("Fitted: ANNUAL = ", fmt(s1$coefficients[1, 1], 3), " + (",
         fmt(s1$coefficients[2, 1], 4), ") * YEAR"),
  paste0("Slope SE = ", fmt(s1$coefficients[2, 2], 4),
         ", t = ", fmt(s1$coefficients[2, 3], 3),
         ", p = ", fmt_p(s1$coefficients[2, 4]), "."),
  paste0("95% CI for slope: [", fmt(ci1[2, 1], 4), ", ", fmt(ci1[2, 2], 4), "] mm per year."),
  paste0("R^2 = ", fmt(s1$r.squared, 4), ", adjusted R^2 = ", fmt(s1$adj.r.squared, 4), "."),
  paste0("Durbin-Watson = ", fmt(dw1, 3), "; residual lag-1 correlation p = ", fmt_p(lag1$p.value), "."),
  paste0("Shapiro-Wilk on residuals p = ", fmt_p(sw_e$p.value),
         "; heteroscedasticity screen (e^2 on fitted) p = ", fmt_p(bp_p), "."),
  paste0("Cook's distance > 4/n flags ", infl_n, " year(s) as influential (not deleted)."),
  "If residuals are autocorrelated, OLS p-values for the slope are not fully reliable;",
  "script 09 reports Mann-Kendall and Sen's slope as rank-based trend inference.",
  "This slope is a descriptive linear summary of 1901-2015, not proof of climate change.",
  "",
  "MODEL 2  Jun-Sep = beta0 + beta1*YEAR + beta2*Jan-Feb + error",
  "Rationale: ask whether winter rainfall is associated with monsoon rainfall after a",
  "linear year term. Jan-Feb is not a subset of Jun-Sep.",
  paste0("R^2 = ", fmt(s2$r.squared, 4), ", adjusted R^2 = ", fmt(s2$adj.r.squared, 4), "."),
  paste0("YEAR coefficient p = ", fmt_p(s2$coefficients["YEAR", 4]),
         "; Jan-Feb coefficient p = ",
         fmt_p(s2$coefficients[grep("Jan", rownames(s2$coefficients))[1], 4]), "."),
  paste0("Durbin-Watson = ", fmt(dw2, 3), ". Approx. VIF (year, Jan-Feb) = ",
         fmt(vif_approx[1], 2), "."),
  "No extra predictors were added merely to raise R^2."
))

save_png("seasonal_correlation_heatmap.png", {
  graphics::par(mar = c(7, 7, 3, 2))
  labs <- c("Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")
  graphics::image(seq_len(4), seq_len(4), t(pear)[4:1, ],
                  col = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(20),
                  zlim = c(-1, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                  main = "Pearson correlation of India seasonal rainfall, 1901-2015")
  graphics::axis(1, at = 1:4, labels = labs, las = 2)
  graphics::axis(2, at = 1:4, labels = rev(labs), las = 1)
  for (i in 1:4) for (j in 1:4) {
    graphics::text(i, 5 - j, sprintf("%.2f", pear[j, i]), cex = 0.9)
  }
})

save_png("scatter_winter_vs_monsoon.png", {
  graphics::plot(india$`Jan-Feb`, india$`Jun-Sep`, pch = 16, col = "#08519C",
                 xlab = "Jan-Feb rainfall (mm)", ylab = "Jun-Sep rainfall (mm)",
                 main = "Winter vs southwest monsoon, India 1901-2015")
  ab <- stats::lm(`Jun-Sep` ~ `Jan-Feb`, data = india)
  graphics::abline(ab, col = "#B30000", lwd = 2)
})

save_png("regression_annual_vs_year.png", {
  graphics::plot(india$YEAR, india$ANNUAL, pch = 16, col = "#08519C",
                 xlab = "Year", ylab = "Annual rainfall (mm)",
                 main = "India annual rainfall vs year with OLS line")
  graphics::abline(m1, col = "#B30000", lwd = 2)
})

save_png("regression_residual_diagnostics.png", {
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
  graphics::plot(f1, e1, pch = 16, col = "#08519C",
                 xlab = "Fitted annual rainfall (mm)", ylab = "Residual (mm)",
                 main = "Residuals vs fitted")
  graphics::abline(h = 0, col = "#B30000", lwd = 2)
  stats::qqnorm(e1, main = "Normal Q-Q of residuals", pch = 16, col = "#08519C")
  stats::qqline(e1, col = "#B30000", lwd = 2)
  graphics::plot(india$YEAR, e1, pch = 16, col = "#08519C",
                 xlab = "Year", ylab = "Residual (mm)", main = "Residuals vs year")
  graphics::abline(h = 0, col = "#B30000", lwd = 2)
  graphics::plot(india$YEAR, cook1, type = "h", col = "#08519C",
                 xlab = "Year", ylab = "Cook's distance",
                 main = "Influence (Cook's distance)")
  graphics::abline(h = 4 / length(cook1), lty = 2, col = "#B30000")
})

cat("07_correlation_regression.R complete.\n")
