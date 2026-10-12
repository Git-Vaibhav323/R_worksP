# Distributional assessment of rainfall. Normality is tested, not assumed.
# KS p-values with estimated parameters are obtained by parametric bootstrap.
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
ann <- finite_vec(india$ANNUAL)
jjas <- finite_vec(india$`Jun-Sep`)

g_ann <- fit_gamma_mle(ann)
n_ann <- fit_normal_mle(ann)
e_ann <- fit_exponential_mle(ann)
g_jjas <- fit_gamma_mle(jjas)

ks_ann_g <- bootstrap_ks(ann, "gamma", B = 400, seed = RNG_SEED)
ks_ann_n <- bootstrap_ks(ann, "norm", B = 400, seed = RNG_SEED + 1L)
ks_ann_e <- bootstrap_ks(ann, "exp", B = 400, seed = RNG_SEED + 2L)
ks_jjas_g <- bootstrap_ks(jjas, "gamma", B = 400, seed = RNG_SEED + 3L)

sw_ann <- stats::shapiro.test(ann)
sw_jjas <- stats::shapiro.test(jjas)

# Ordinary KS is reported only as a diagnostic statistic; p-value is NOT used.
ks_ordinary_note <- stats::ks.test(ann, "pgamma", shape = g_ann$shape, rate = g_ann$rate)

fit_tbl <- rbind(
  data.frame(series = "India annual", distribution = "gamma",
             param1_name = "shape", param1 = g_ann$shape,
             param2_name = "rate", param2 = g_ann$rate,
             loglik = g_ann$loglik, aic = g_ann$aic,
             ks_stat = ks_ann_g$statistic, ks_p_bootstrap = ks_ann_g$p_bootstrap,
             shapiro_W = sw_ann$statistic, shapiro_p = sw_ann$p.value,
             stringsAsFactors = FALSE),
  data.frame(series = "India annual", distribution = "normal",
             param1_name = "mean", param1 = n_ann$mean,
             param2_name = "sd_mle", param2 = n_ann$sd_mle,
             loglik = n_ann$loglik, aic = n_ann$aic,
             ks_stat = ks_ann_n$statistic, ks_p_bootstrap = ks_ann_n$p_bootstrap,
             shapiro_W = sw_ann$statistic, shapiro_p = sw_ann$p.value,
             stringsAsFactors = FALSE),
  data.frame(series = "India annual", distribution = "exponential",
             param1_name = "rate", param1 = e_ann$rate,
             param2_name = "mean", param2 = e_ann$mean,
             loglik = e_ann$loglik, aic = e_ann$aic,
             ks_stat = ks_ann_e$statistic, ks_p_bootstrap = ks_ann_e$p_bootstrap,
             shapiro_W = sw_ann$statistic, shapiro_p = sw_ann$p.value,
             stringsAsFactors = FALSE),
  data.frame(series = "India Jun-Sep", distribution = "gamma",
             param1_name = "shape", param1 = g_jjas$shape,
             param2_name = "rate", param2 = g_jjas$rate,
             loglik = g_jjas$loglik, aic = g_jjas$aic,
             ks_stat = ks_jjas_g$statistic, ks_p_bootstrap = ks_jjas_g$p_bootstrap,
             shapiro_W = sw_jjas$statistic, shapiro_p = sw_jjas$p.value,
             stringsAsFactors = FALSE)
)
write.csv(fit_tbl, file.path("results", "tables", "distribution_fit_summary.csv"), row.names = FALSE)

# Two-part description for a dry-region winter month (zeros are real)
subdiv <- dat$subdiv
if (is.character(subdiv$months_complete)) {
  wr <- subdiv[subdiv$SUBDIVISION == "West Rajasthan" & subdiv$months_complete %in% c("TRUE", "true"), ]
} else {
  wr <- subdiv[subdiv$SUBDIVISION == "West Rajasthan" & subdiv$months_complete, ]
}
dec <- finite_vec(wr$DEC)
n_dec <- length(dec)
n_zero <- sum(dec == 0)
pos <- dec[dec > 0]
g_dec <- if (length(pos) >= 5L) fit_gamma_mle(pos) else NULL
zero_fit <- data.frame(
  series = "West Rajasthan December (complete years)",
  n = n_dec, n_zero = n_zero, p_zero = n_zero / n_dec,
  n_positive = length(pos),
  gamma_shape_positive = if (!is.null(g_dec)) g_dec$shape else NA_real_,
  gamma_rate_positive = if (!is.null(g_dec)) g_dec$rate else NA_real_,
  note = "Hurdle: Bernoulli(zero) plus gamma on positive amounts. Exponential is not used for the zero mass.",
  stringsAsFactors = FALSE
)
write.csv(zero_fit, file.path("results", "tables", "zero_inflated_month_example.csv"), row.names = FALSE)

write_txt(file.path("results", "statistical_outputs", "distribution_fit_notes.txt"), c(
  "PROBABILITY DISTRIBUTIONS",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Bootstrap replications B =", ks_ann_g$B, "; seed", RNG_SEED),
  "",
  "WHY NOT ASSUME NORMALITY",
  "Annual rainfall is bounded below by 0 and typically right-skewed. A normal model",
  "assigns probability to negative rainfall. Shapiro-Wilk tests the normal model directly.",
  paste0("India annual Shapiro-Wilk W = ", fmt(as.numeric(sw_ann$statistic), 4),
         ", p = ", fmt_p(sw_ann$p.value), "."),
  "",
  "GAMMA",
  "A gamma random variable is positive and can be right-skewed. Parameters were estimated",
  "by maximum likelihood (nlm), started at method-of-moments values.",
  paste0("Annual gamma MLE: shape = ", fmt(g_ann$shape, 3), ", rate = ", fmt(g_ann$rate, 5),
         ", scale = ", fmt(g_ann$scale, 2), " mm; AIC = ", fmt(g_ann$aic, 1), "."),
  paste0("Bootstrap KS (parameters re-estimated in each replicate): D = ",
         fmt(ks_ann_g$statistic, 4), ", p = ", fmt_p(ks_ann_g$p_bootstrap), "."),
  "The ordinary ks.test p-value with plugged-in MLE is NOT valid (Lilliefors problem) and is not used.",
  paste0("For reference only, unadjusted ks.test D = ", fmt(as.numeric(ks_ordinary_note$statistic), 4), "."),
  "",
  "EXPONENTIAL",
  "The exponential law is a gamma with shape = 1. It is a waiting-time model with",
  "memorylessness and the restriction SD = mean. All-India annual rainfall has mean near",
  paste0(fmt(mean(ann), 1), " mm and SD near ", fmt(stats::sd(ann), 1),
         " mm, so SD << mean; exponential is scientifically implausible here."),
  "It is fitted only to demonstrate a poor reference model, not as a candidate climate law.",
  paste0("Exponential AIC = ", fmt(e_ann$aic, 1), " versus gamma AIC = ", fmt(g_ann$aic, 1),
         " and normal AIC = ", fmt(n_ann$aic, 1), "."),
  paste0("Bootstrap KS for exponential: D = ", fmt(ks_ann_e$statistic, 4),
         ", p = ", fmt_p(ks_ann_e$p_bootstrap), "."),
  "",
  "ZEROS",
  "National monthly and annual series have no zeros. Subdivision dry-season months do.",
  "A single gamma or exponential cannot place a point mass at zero. The West Rajasthan",
  "December example uses an empirical P(zero) plus gamma on positive millimetres.",
  paste0("West Rajasthan DEC: n = ", n_dec, ", zeros = ", n_zero,
         ", P(zero) = ", fmt(n_zero / n_dec, 3), "."),
  "",
  "LIMITATIONS",
  "Years are temporally dependent, so iid goodness-of-fit p-values are approximate.",
  "Bootstrap B = 400 is adequate for reporting to two decimals, not for tiny p-values."
))

# Plots
save_png("india_annual_histogram_density.png", {
  graphics::hist(ann, breaks = 15, freq = FALSE, col = "#C6DBEF", border = "white",
                 xlab = "Annual rainfall (mm)", ylab = "Density",
                 main = "India annual rainfall: histogram and kernel density, 1901-2015")
  graphics::lines(stats::density(ann), lwd = 2, col = "#08306B")
  graphics::abline(v = mean(ann), col = "#B30000", lwd = 2, lty = 2)
  graphics::abline(v = stats::median(ann), col = "#006D2C", lwd = 2, lty = 3)
  graphics::legend("topright", c("Kernel density", "Mean", "Median"),
                   col = c("#08306B", "#B30000", "#006D2C"), lwd = 2, lty = c(1, 2, 3),
                   bty = "n", cex = 0.8)
})

save_png("india_annual_ecdf.png", {
  graphics::plot(stats::ecdf(ann), xlab = "Annual rainfall (mm)", ylab = "F(x)",
                 main = "Empirical CDF of India annual rainfall, 1901-2015",
                 col = "#08306B", lwd = 2)
  xx <- seq(min(ann), max(ann), length.out = 200)
  graphics::lines(xx, stats::pgamma(xx, shape = g_ann$shape, rate = g_ann$rate),
                  col = "#B30000", lwd = 2, lty = 2)
  graphics::legend("bottomright", c("Empirical CDF", "Fitted gamma CDF"),
                   col = c("#08306B", "#B30000"), lwd = 2, lty = c(1, 2), bty = "n")
})

save_png("india_annual_qq_normal.png", {
  stats::qqnorm(ann, main = "Normal Q-Q plot: India annual rainfall (mm)",
                xlab = "Theoretical normal quantiles", ylab = "Observed annual rainfall (mm)",
                pch = 16, col = "#08519C")
  stats::qqline(ann, col = "#B30000", lwd = 2)
})

save_png("india_annual_qq_gamma.png", {
  p <- stats::ppoints(length(ann))
  th <- stats::qgamma(p, shape = g_ann$shape, rate = g_ann$rate)
  graphics::plot(sort(th), sort(ann), pch = 16, col = "#08519C",
                 xlab = "Theoretical gamma quantiles (mm)",
                 ylab = "Observed annual rainfall (mm)",
                 main = "Gamma Q-Q plot: India annual rainfall")
  graphics::abline(0, 1, col = "#B30000", lwd = 2)
})

save_png("india_annual_distribution_overlay.png", {
  graphics::hist(ann, breaks = 15, freq = FALSE, col = "#F0F0F0", border = "grey70",
                 xlab = "Annual rainfall (mm)", ylab = "Density",
                 main = "Fitted densities vs India annual rainfall, 1901-2015")
  xx <- seq(min(ann) - 50, max(ann) + 50, length.out = 300)
  graphics::lines(xx, stats::dgamma(xx, g_ann$shape, g_ann$rate), col = "#0072B2", lwd = 2)
  graphics::lines(xx, stats::dnorm(xx, n_ann$mean, n_ann$sd_mle), col = "#D55E00", lwd = 2, lty = 2)
  graphics::lines(xx, stats::dexp(xx, e_ann$rate), col = "#009E73", lwd = 2, lty = 3)
  graphics::legend("topright", c("Gamma MLE", "Normal MLE", "Exponential MLE"),
                   col = c("#0072B2", "#D55E00", "#009E73"), lwd = 2, lty = c(1, 2, 3),
                   bty = "n", cex = 0.8)
})

cat("04_probability_distributions.R complete.\n")
