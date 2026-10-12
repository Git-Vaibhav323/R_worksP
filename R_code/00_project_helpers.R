# Shared utilities for Phase 2 analysis. Base R only.
# Sourced by 03–09. Sets the project root when run via Rscript.

ensure_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    d <- dirname(script_path)
    if (basename(d) == "R_code") setwd(dirname(d)) else setwd(d)
  } else if (basename(getwd()) == "R_code") {
    setwd("..")
  }
  if (!dir.exists("data") || !dir.exists("results")) {
    stop("Project root not found. Set the working directory to the folder that contains data/ and R_code/.")
  }
  invisible(getwd())
}

ensure_project_root()

MONTH_COLS <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
MONTH_LABELS <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# IMD-style monsoon categories relative to a long-period average (LPA).
# Here LPA is the sample mean of the series under study, documented as such.
LPA_DEFICIENT <- 0.90
LPA_EXCESS <- 1.10

ALPHA <- 0.05
RNG_SEED <- 19012015

finite_vec <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

skewness_g1 <- function(x) {
  x <- finite_vec(x)
  n <- length(x)
  if (n < 3L) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  n / ((n - 1) * (n - 2)) * sum(((x - m) / s)^3)
}

kurtosis_g2 <- function(x) {
  # Excess kurtosis (normal = 0), unbiased estimator.
  x <- finite_vec(x)
  n <- length(x)
  if (n < 4L) return(NA_real_)
  m <- mean(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NA_real_)
  z <- (x - m) / s
  a <- n * (n + 1) / ((n - 1) * (n - 2) * (n - 3)) * sum(z^4)
  b <- 3 * (n - 1)^2 / ((n - 2) * (n - 3))
  a - b
}

sample_mode <- function(x, digits = 1) {
  x <- finite_vec(x)
  if (!length(x)) return(NA_real_)
  xr <- round(x, digits)
  tab <- sort(table(xr), decreasing = TRUE)
  if (length(tab) == 0L) return(NA_real_)
  if (length(tab) > 1L && tab[1] == tab[2]) return(NA_real_)
  as.numeric(names(tab)[1])
}

desc_stats <- function(x, series_name, unit = "mm") {
  x_all <- as.numeric(x)
  n_miss <- sum(!is.finite(x_all))
  x <- finite_vec(x_all)
  n <- length(x)
  if (!n) {
    return(data.frame(
      series = series_name, unit = unit, n = 0L, n_missing = n_miss,
      mean = NA_real_, median = NA_real_, mode_rounded = NA_real_,
      min = NA_real_, max = NA_real_, range = NA_real_,
      variance = NA_real_, sd = NA_real_, se = NA_real_,
      q1 = NA_real_, q3 = NA_real_, iqr = NA_real_,
      p5 = NA_real_, p10 = NA_real_, p90 = NA_real_, p95 = NA_real_,
      cv_pct = NA_real_, skewness = NA_real_, excess_kurtosis = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  qs <- stats::quantile(x, probs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
                        names = FALSE, type = 7)
  m <- mean(x)
  s <- stats::sd(x)
  data.frame(
    series = series_name,
    unit = unit,
    n = n,
    n_missing = n_miss,
    mean = m,
    median = qs[4],
    mode_rounded = sample_mode(x, 1),
    min = min(x),
    max = max(x),
    range = max(x) - min(x),
    variance = stats::var(x),
    sd = s,
    se = s / sqrt(n),
    q1 = qs[3],
    q3 = qs[5],
    iqr = qs[5] - qs[3],
    p5 = qs[1],
    p10 = qs[2],
    p90 = qs[6],
    p95 = qs[7],
    cv_pct = if (m != 0) 100 * s / m else NA_real_,
    skewness = skewness_g1(x),
    excess_kurtosis = kurtosis_g2(x),
    stringsAsFactors = FALSE
  )
}

save_png <- function(filename, expr, width = 2400, height = 1600, res = 300) {
  path <- file.path("graphs", filename)
  grDevices::png(path, width = width, height = height, res = res)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3.2, 1.2), cex.main = 0.95, cex.lab = 0.9, cex.axis = 0.8)
  force(expr)
  invisible(path)
}

write_txt <- function(path, lines) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat(paste(lines, collapse = "\n"), "\n", file = con, sep = "")
}

append_txt <- function(path, lines) {
  cat(paste(lines, collapse = "\n"), "\n", file = path, append = TRUE, sep = "")
}

load_cleaned <- function() {
  india <- read.csv(
    file.path("results", "cleaned_data", "india_national_rainfall_cleaned.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  subdiv <- read.csv(
    file.path("results", "cleaned_data", "imd_subdivision_rainfall_cleaned.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  india$YEAR <- as.integer(india$YEAR)
  subdiv$YEAR <- as.integer(subdiv$YEAR)
  subdiv$months_complete <- as.logical(subdiv$months_complete)
  list(india = india, subdiv = subdiv)
}

monsoon_category <- function(x, lpa, lo = LPA_DEFICIENT, hi = LPA_EXCESS) {
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x) & is.finite(lpa)
  out[ok & x < lo * lpa] <- "deficient"
  out[ok & x >= lo * lpa & x <= hi * lpa] <- "normal"
  out[ok & x > hi * lpa] <- "excess"
  out
}

# Gamma MLE (shape, rate). Moments start; nlm on negative log-likelihood.
fit_gamma_mle <- function(x) {
  x <- finite_vec(x)
  x <- x[x > 0]
  n <- length(x)
  stopifnot(n >= 5L)
  m <- mean(x)
  v <- stats::var(x)
  shape0 <- max(m^2 / v, 0.1)
  rate0 <- max(shape0 / m, 1e-6)
  nll <- function(par) {
    sh <- par[1]; ra <- par[2]
    if (sh <= 0 || ra <= 0) return(1e12)
    -sum(stats::dgamma(x, shape = sh, rate = ra, log = TRUE))
  }
  fit <- stats::nlm(nll, c(shape0, rate0), steptol = 1e-8, gradtol = 1e-8, iterlim = 200)
  sh <- fit$estimate[1]
  ra <- fit$estimate[2]
  ll <- -nll(c(sh, ra))
  list(
    n = n, mean = m, var = v,
    shape_mom = shape0, rate_mom = rate0,
    shape = sh, rate = ra, scale = 1 / ra,
    loglik = ll, aic = -2 * ll + 4, nlm_code = fit$code
  )
}

fit_exponential_mle <- function(x) {
  x <- finite_vec(x)
  x <- x[x > 0]
  n <- length(x)
  rate <- 1 / mean(x)
  ll <- sum(stats::dexp(x, rate = rate, log = TRUE))
  list(n = n, rate = rate, mean = 1 / rate, loglik = ll, aic = -2 * ll + 2)
}

fit_normal_mle <- function(x) {
  x <- finite_vec(x)
  n <- length(x)
  mu <- mean(x)
  # MLE sd uses n, not n-1
  sig <- sqrt(mean((x - mu)^2))
  ll <- sum(stats::dnorm(x, mu, sig, log = TRUE))
  list(n = n, mean = mu, sd_mle = sig, sd_unbiased = stats::sd(x),
       loglik = ll, aic = -2 * ll + 4)
}

# Parametric bootstrap KS (Lilliefors-type) for a named distribution.
bootstrap_ks <- function(x, dist = c("gamma", "norm", "exp"), B = 400, seed = RNG_SEED) {
  dist <- match.arg(dist)
  x <- finite_vec(x)
  if (dist != "norm") x <- x[x > 0]
  n <- length(x)
  set.seed(seed)
  ks_stat <- function(sample, d) {
    if (d == "gamma") {
      f <- fit_gamma_mle(sample)
      stats::ks.test(sample, "pgamma", shape = f$shape, rate = f$rate)$statistic
    } else if (d == "norm") {
      f <- fit_normal_mle(sample)
      stats::ks.test(sample, "pnorm", mean = f$mean, sd = f$sd_mle)$statistic
    } else {
      f <- fit_exponential_mle(sample)
      stats::ks.test(sample, "pexp", rate = f$rate)$statistic
    }
  }
  obs <- as.numeric(ks_stat(x, dist))
  sim <- numeric(B)
  for (b in seq_len(B)) {
    if (dist == "gamma") {
      f <- fit_gamma_mle(x)
      y <- stats::rgamma(n, shape = f$shape, rate = f$rate)
    } else if (dist == "norm") {
      f <- fit_normal_mle(x)
      y <- stats::rnorm(n, f$mean, f$sd_mle)
    } else {
      f <- fit_exponential_mle(x)
      y <- stats::rexp(n, rate = f$rate)
    }
    sim[b] <- tryCatch(as.numeric(ks_stat(y, dist)), error = function(e) NA_real_)
  }
  sim <- sim[is.finite(sim)]
  p <- mean(sim >= obs)
  list(statistic = obs, p_bootstrap = p, B = length(sim), n = n, dist = dist)
}

mann_kendall <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  x <- x[ok]
  n <- length(x)
  if (n < 8L) stop("Mann-Kendall needs a longer series.")
  S <- 0
  for (i in seq_len(n - 1L)) {
    S <- S + sum(sign(x[(i + 1L):n] - x[i]))
  }
  # Tie correction
  tab <- table(x)
  ties <- as.numeric(tab[tab > 1])
  varS <- n * (n - 1) * (2 * n + 5)
  if (length(ties)) {
    varS <- varS - sum(ties * (ties - 1) * (2 * ties + 5))
  }
  varS <- varS / 18
  if (S > 0) Z <- (S - 1) / sqrt(varS)
  else if (S < 0) Z <- (S + 1) / sqrt(varS)
  else Z <- 0
  p <- 2 * stats::pnorm(-abs(Z))
  list(n = n, S = S, varS = varS, Z = Z, p_value = p,
       tau = S / (n * (n - 1) / 2))
}

sens_slope <- function(x, t = NULL) {
  x <- as.numeric(x)
  if (is.null(t)) t <- seq_along(x)
  ok <- is.finite(x) & is.finite(t)
  x <- x[ok]; t <- t[ok]
  n <- length(x)
  slopes <- numeric(n * (n - 1) / 2)
  k <- 0L
  for (i in seq_len(n - 1L)) {
    dt <- t[(i + 1L):n] - t[i]
    dx <- x[(i + 1L):n] - x[i]
    keep <- dt != 0
    m <- sum(keep)
    if (!m) next
    slopes[(k + 1):(k + m)] <- dx[keep] / dt[keep]
    k <- k + m
  }
  slopes <- slopes[seq_len(k)]
  med <- stats::median(slopes)
  list(n = n, n_pairs = length(slopes), slope = med,
       intercept = stats::median(x - med * t))
}

durbin_watson <- function(e) {
  e <- finite_vec(e)
  sum(diff(e)^2) / sum(e^2)
}

fmt <- function(x, d = 3) {
  ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "NA")
}

fmt_p <- function(p) {
  vapply(p, function(pi) {
    if (!is.finite(pi)) return("NA")
    if (pi < 0.0001) return("< 0.0001")
    sprintf("%.4f", pi)
  }, character(1))
}
