# Time-series description and trend inference. Not a climate-change attribution study.
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
india <- india[order(india$YEAR), ]
ann <- india$ANNUAL
jjas <- india$`Jun-Sep`
year <- india$YEAR
n <- length(ann)

ma <- function(x, k = 11) {
  stats::filter(x, rep(1 / k, k), sides = 2)
}
ann_ma11 <- as.numeric(ma(ann, 11))
jjas_ma11 <- as.numeric(ma(jjas, 11))

mk_ann <- mann_kendall(ann)
mk_jjas <- mann_kendall(jjas)
sen_ann <- sens_slope(ann, year)
sen_jjas <- sens_slope(jjas, year)

acf_ann <- stats::acf(ann, lag.max = 15, plot = FALSE)
pacf_ann <- stats::pacf(ann, lag.max = 15, plot = FALSE)
acf_e <- stats::acf(stats::residuals(stats::lm(ann ~ year)), lag.max = 15, plot = FALSE)

m_ols <- stats::lm(ann ~ year)
ols_s <- summary(m_ols)
dw <- durbin_watson(stats::residuals(m_ols))

# Approximate 95% bands for ACF under white noise: +/- 1.96/sqrt(n)
acf_tbl <- data.frame(
  lag = as.numeric(acf_ann$lag),
  acf_annual = as.numeric(acf_ann$acf),
  acf_ols_residual = as.numeric(acf_e$acf),
  white_noise_bound = 1.96 / sqrt(n)
)
write.csv(acf_tbl, file.path("results", "tables", "autocorrelation_annual.csv"), row.names = FALSE)

trend_tbl <- data.frame(
  series = c("India annual", "India Jun-Sep"),
  n = c(mk_ann$n, mk_jjas$n),
  ols_slope_mm_per_year = c(unname(coef(m_ols)[2]), unname(coef(stats::lm(jjas ~ year))[2])),
  ols_slope_p = c(ols_s$coefficients[2, 4], summary(stats::lm(jjas ~ year))$coefficients[2, 4]),
  pearson_r_with_year = c(stats::cor(year, ann), stats::cor(year, jjas)),
  mann_kendall_S = c(mk_ann$S, mk_jjas$S),
  mann_kendall_Z = c(mk_ann$Z, mk_jjas$Z),
  mann_kendall_tau = c(mk_ann$tau, mk_jjas$tau),
  mann_kendall_p = c(mk_ann$p_value, mk_jjas$p_value),
  sen_slope_mm_per_year = c(sen_ann$slope, sen_jjas$slope),
  sen_intercept = c(sen_ann$intercept, sen_jjas$intercept),
  durbin_watson_ols = c(dw, durbin_watson(stats::residuals(stats::lm(jjas ~ year)))),
  stringsAsFactors = FALSE
)
write.csv(trend_tbl, file.path("results", "tables", "trend_tests_national.csv"), row.names = FALSE)

ma_tbl <- data.frame(YEAR = year, ANNUAL = ann, annual_ma11 = ann_ma11,
                     JJAS = jjas, jjas_ma11 = jjas_ma11)
write.csv(ma_tbl, file.path("results", "tables", "annual_moving_average.csv"), row.names = FALSE)

# Regional MK for selected complete series (same six as ANOVA)
subdiv <- dat$subdiv
if (is.character(subdiv$months_complete)) {
  sub_ok <- subdiv[subdiv$months_complete %in% c("TRUE", "true"), ]
} else {
  sub_ok <- subdiv[subdiv$months_complete == TRUE, ]
}
sub_ok <- sub_ok[sub_ok$YEAR >= 1901 & sub_ok$YEAR <= 2015, ]
selected <- c("Kerala", "Coastal Karnataka", "Gangetic West Bengal",
              "Tamil Nadu", "Punjab", "West Rajasthan")
reg_mk <- lapply(selected, function(rg) {
  d <- sub_ok[sub_ok$SUBDIVISION == rg, ]
  d <- d[order(d$YEAR), ]
  mk <- mann_kendall(d$ANNUAL)
  se <- sens_slope(d$ANNUAL, d$YEAR)
  data.frame(
    SUBDIVISION = rg, n = mk$n,
    mk_tau = mk$tau, mk_Z = mk$Z, mk_p = mk$p_value,
    sen_slope_mm_per_year = se$slope,
    stringsAsFactors = FALSE
  )
})
write.csv(do.call(rbind, reg_mk),
          file.path("results", "tables", "regional_mann_kendall.csv"), row.names = FALSE)

write_txt(file.path("results", "statistical_outputs", "time_series_trend_notes.txt"), c(
  "TIME-SERIES AND TREND ANALYSIS",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "THIS IS NOT A CLIMATE-CHANGE ATTRIBUTION STUDY.",
  "A slope or Mann-Kendall p-value describes monotonic association with time in",
  "1901-2015. It does not identify greenhouse gases, aerosols, land-use, or natural",
  "decadal variability as causes, and it does not forecast the next century.",
  "",
  "DESCRIPTIVE TOOLS",
  "The 11-year moving average is a smoother (centred). End years are missing because",
  "the window is incomplete. It is not a statistical test.",
  "",
  "OLS TREND vs RANK-BASED TREND",
  paste0("Annual OLS slope = ", fmt(unname(coef(m_ols)[2]), 4), " mm/year, p = ",
         fmt_p(ols_s$coefficients[2, 4]), ", R^2 = ", fmt(ols_s$r.squared, 4), "."),
  paste0("Durbin-Watson = ", fmt(dw, 3), ". Values near 2 suggest little lag-1 residual correlation;",
         " values << 2 suggest positive autocorrelation and anti-conservative OLS p-values."),
  paste0("Mann-Kendall on annual totals: S = ", mk_ann$S, ", tau = ", fmt(mk_ann$tau, 3),
         ", Z = ", fmt(mk_ann$Z, 3), ", p = ", fmt_p(mk_ann$p_value), "."),
  paste0("Sen's slope (median pairwise slope) = ", fmt(sen_ann$slope, 4), " mm/year."),
  paste0("Jun-Sep Mann-Kendall p = ", fmt_p(mk_jjas$p_value),
         ", Sen slope = ", fmt(sen_jjas$slope, 4), " mm/year."),
  "Mann-Kendall assumes independent observations under H0 (or uses the simple variance",
  "formula implemented here). Strong autocorrelation would require a modified variance;",
  "inspect the ACF before treating a borderline p-value as decisive.",
  "",
  "AUTOCORRELATION",
  paste0("Lag-1 sample ACF of annual rainfall = ", fmt(as.numeric(acf_ann$acf[2]), 3),
         " (white-noise 95% bound about ", fmt(1.96 / sqrt(n), 3), ")."),
  "If |ACF| stays inside the bound, serial dependence is not dominating the series.",
  "",
  "DECISION LANGUAGE",
  paste("At alpha = 0.05, annual Mann-Kendall:",
        ifelse(mk_ann$p_value < ALPHA,
               "reject no-trend (monotonic association with time detected).",
               "do not reject no-trend; absence of a significant MK statistic is not proof of stationarity.")),
  "Distinguish this from the OLS line, which always exists as a descriptive fit."
))

save_png("india_annual_timeseries.png", {
  graphics::plot(year, ann, type = "l", col = "#4D4D4D", lwd = 1,
                 xlab = "Year", ylab = "Annual rainfall (mm)",
                 main = "India area-weighted annual rainfall, 1901-2015")
  graphics::points(year, ann, pch = 16, cex = 0.4, col = "#08519C")
  graphics::lines(year, ann_ma11, col = "#B30000", lwd = 2)
  graphics::abline(m_ols, col = "#1B9E77", lwd = 2, lty = 2)
  graphics::legend("topright",
                   c("Annual series", "11-year moving average", "OLS line"),
                   col = c("#4D4D4D", "#B30000", "#1B9E77"), lwd = c(1, 2, 2),
                   lty = c(1, 1, 2), bty = "n", cex = 0.8)
})

save_png("india_jjas_timeseries.png", {
  graphics::plot(year, jjas, type = "l", col = "#4D4D4D",
                 xlab = "Year", ylab = "Jun-Sep rainfall (mm)",
                 main = "India southwest monsoon rainfall, 1901-2015")
  graphics::lines(year, jjas_ma11, col = "#B30000", lwd = 2)
  graphics::legend("topright", c("JJAS", "11-year MA"), col = c("#4D4D4D", "#B30000"),
                   lwd = c(1, 2), bty = "n")
})

save_png("india_annual_acf.png", {
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
  stats::acf(ann, lag.max = 15, main = "ACF: India annual rainfall")
  stats::pacf(ann, lag.max = 15, main = "PACF: India annual rainfall")
})

cat("09_time_series_analysis.R complete.\n")
