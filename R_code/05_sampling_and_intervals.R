# Sampling distributions, CLT demonstration, and confidence intervals.
# Simulated draws are teaching devices, not extra observed rainfall years.
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
ann <- finite_vec(india$ANNUAL)
jjas <- finite_vec(india$`Jun-Sep`)
n <- length(ann)
xbar <- mean(ann)
s <- stats::sd(ann)
se <- s / sqrt(n)

mean_ci <- function(x, conf = 0.95) {
  x <- finite_vec(x)
  n <- length(x)
  m <- mean(x)
  se <- stats::sd(x) / sqrt(n)
  tcrit <- stats::qt(1 - (1 - conf) / 2, df = n - 1)
  c(mean = m, se = se, df = n - 1, conf = conf,
    lower = m - tcrit * se, upper = m + tcrit * se)
}

prop_ci <- function(successes, n, conf = 0.95) {
  p <- successes / n
  z <- stats::qnorm(1 - (1 - conf) / 2)
  se <- sqrt(p * (1 - p) / n)
  c(p = p, n = n, successes = successes, conf = conf, se = se,
    lower = max(0, p - z * se), upper = min(1, p + z * se))
}

levels <- c(0.90, 0.95, 0.99)
ci_mean <- do.call(rbind, lapply(levels, function(lv) {
  v <- mean_ci(ann, lv)
  data.frame(series = "India annual rainfall",
             parameter = "mean (mm)",
             confidence_level = lv, n = n, estimate = v["mean"], se = v["se"],
             df = v["df"], lower = v["lower"], upper = v["upper"],
             method = "Student t, iid assumption",
             stringsAsFactors = FALSE)
}))
ci_mean_jjas <- do.call(rbind, lapply(levels, function(lv) {
  v <- mean_ci(jjas, lv)
  nj <- length(jjas)
  data.frame(series = "India Jun-Sep rainfall",
             parameter = "mean (mm)",
             confidence_level = lv, n = nj, estimate = v["mean"], se = v["se"],
             df = v["df"], lower = v["lower"], upper = v["upper"],
             method = "Student t, iid assumption",
             stringsAsFactors = FALSE)
}))

lpa <- mean(jjas)
n_def <- sum(jjas < LPA_DEFICIENT * lpa)
ci_prop <- do.call(rbind, lapply(levels, function(lv) {
  v <- prop_ci(n_def, n, lv)
  data.frame(series = "India Jun-Sep deficient years",
             parameter = "proportion deficient (<90% of series mean)",
             confidence_level = lv, n = n, estimate = v["p"], se = v["se"],
             df = NA_real_, lower = v["lower"], upper = v["upper"],
             method = "Wald (normal) interval for a binomial proportion",
             stringsAsFactors = FALSE)
}))
ci_tbl <- rbind(ci_mean, ci_mean_jjas, ci_prop)
write.csv(ci_tbl, file.path("results", "tables", "confidence_intervals.csv"), row.names = FALSE)

# CLT simulation: resampling from the empirical annual distribution
set.seed(RNG_SEED)
n_rep <- 4000
sizes <- c(5L, 15L, 40L)
sim_rows <- list()
sim_means <- list()
for (ni in sizes) {
  mns <- replicate(n_rep, mean(sample(ann, size = ni, replace = TRUE)))
  sim_means[[as.character(ni)]] <- mns
  sim_rows[[length(sim_rows) + 1L]] <- data.frame(
    sample_size = ni,
    n_replications = n_rep,
    mean_of_sample_means = mean(mns),
    sd_of_sample_means = stats::sd(mns),
    theoretical_se = s / sqrt(ni),
    skewness_of_means = skewness_g1(mns),
    shapiro_p_500 = stats::shapiro.test(mns[seq_len(500)])$p.value,
    stringsAsFactors = FALSE
  )
}
sim_tbl <- do.call(rbind, sim_rows)
write.csv(sim_tbl, file.path("results", "tables", "clt_simulation_summary.csv"), row.names = FALSE)

# Sampling distribution vs parent
save_png("clt_sampling_distributions.png", {
  graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3, 1))
  graphics::hist(ann, breaks = 12, col = "#C6DBEF", main = "Parent: annual rainfall (n = 115 years)",
                 xlab = "mm", ylab = "Frequency")
  for (ni in sizes) {
    graphics::hist(sim_means[[as.character(ni)]], breaks = 30, col = "#9ECAE1",
                   main = paste0("Means of resamples, n = ", ni),
                   xlab = "Sample mean (mm)", ylab = "Frequency")
    graphics::abline(v = xbar, col = "#B30000", lwd = 2)
  }
})

save_png("confidence_intervals_mean_compare.png", {
  sub <- ci_mean
  ys <- seq_len(nrow(sub))
  graphics::plot(sub$estimate, ys, xlim = range(c(sub$lower, sub$upper)),
                 yaxt = "n", xlab = "Annual mean rainfall (mm)", ylab = "",
                 main = "t intervals for India mean annual rainfall",
                 pch = 16, col = "#08306B")
  graphics::axis(2, at = ys, labels = paste0(100 * sub$confidence_level, "%"), las = 1)
  graphics::arrows(sub$lower, ys, sub$upper, ys, angle = 90, code = 3, length = 0.05, lwd = 2,
                   col = "#08519C")
  graphics::abline(v = xbar, lty = 2, col = "grey40")
})

write_txt(file.path("results", "statistical_outputs", "sampling_and_intervals_notes.txt"), c(
  "SAMPLING, CENTRAL LIMIT THEOREM, AND CONFIDENCE INTERVALS",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Random seed:", RNG_SEED),
  "",
  "TEACHING SIMULATION VERSUS DATA",
  "The CLT panels resample WITH REPLACEMENT from the 115 observed annual totals.",
  "They are not additional meteorological years and must not be averaged into the data.",
  "",
  "CLT DESIGN",
  paste0("Parent mean = ", fmt(xbar, 2), " mm, parent SD = ", fmt(s, 2), " mm."),
  paste0("Replications = ", n_rep, " for sample sizes ", paste(sizes, collapse = ", "), "."),
  "As n grows, the histogram of sample means concentrates at the parent mean and",
  "becomes more symmetric. Theoretical SE = parent_SD / sqrt(n).",
  paste0("For n = 40, simulated SD of means = ", fmt(sim_tbl$sd_of_sample_means[sim_tbl$sample_size == 40], 3),
         " mm; theoretical SE = ", fmt(sim_tbl$theoretical_se[sim_tbl$sample_size == 40], 3), " mm."),
  "",
  "CONFIDENCE INTERVALS",
  "Mean intervals use Student's t with df = n-1. They assume iid normal (or approximately",
  "normal means by the CLT). Serial correlation in annual rainfall can make the true",
  "uncertainty larger than these intervals; treat them as iid-benchmark intervals.",
  paste0("95% CI for mean annual rainfall: ",
         fmt(ci_mean$lower[ci_mean$confidence_level == 0.95], 2), " to ",
         fmt(ci_mean$upper[ci_mean$confidence_level == 0.95], 2), " mm."),
  "Wald intervals for the deficient-monsoon proportion assume iid Bernoulli trials;",
  "years are not independent experiments. Higher confidence levels produce wider intervals.",
  "",
  "INFERENCE FROM A HISTORICAL TIME SERIES",
  "The 115 years are one realised path of a climate process, not a simple random sample",
  "from a super-population of years. Estimates describe this record. Probability statements",
  "about a 'random year' require extra assumptions (stationarity, weak dependence)."
))

cat("05_sampling_and_intervals.R complete.\n")
