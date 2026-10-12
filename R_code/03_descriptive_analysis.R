# Descriptive statistics and empirical probabilities for IMD rainfall (mm).
source("R_code/00_project_helpers.R")

dat <- load_cleaned()
india <- dat$india
subdiv <- dat$subdiv
sub_ok <- subdiv[isTRUE(subdiv$months_complete) | subdiv$months_complete %in% TRUE, ]
# read.csv may leave months_complete as TRUE/FALSE character already converted
if (is.character(subdiv$months_complete)) {
  sub_ok <- subdiv[subdiv$months_complete %in% c("TRUE", "True", "true"), ]
} else {
  sub_ok <- subdiv[subdiv$months_complete, ]
}

rows <- list()
push <- function(df) rows[[length(rows) + 1L]] <<- df

push(desc_stats(india$ANNUAL, "India annual (1901-2015)"))
push(desc_stats(india$`Jan-Feb`, "India Jan-Feb winter"))
push(desc_stats(india$`Mar-May`, "India Mar-May pre-monsoon"))
push(desc_stats(india$`Jun-Sep`, "India Jun-Sep southwest monsoon"))
push(desc_stats(india$`Oct-Dec`, "India Oct-Dec post-monsoon"))
for (m in MONTH_COLS) {
  push(desc_stats(india[[m]], paste("India", m)))
}

# Regional annual, complete years, common national overlap 1901-2015
sub_common <- sub_ok[sub_ok$YEAR >= 1901 & sub_ok$YEAR <= 2015, ]
regs <- sort(unique(sub_common$SUBDIVISION))
for (rg in regs) {
  d <- sub_common[sub_common$SUBDIVISION == rg, ]
  push(desc_stats(d$ANNUAL, paste0("Subdivision annual: ", rg, " (complete years 1901-2015)")))
  push(desc_stats(d$JJAS, paste0("Subdivision JJAS: ", rg)))
}

desc_tbl <- do.call(rbind, rows)
write.csv(desc_tbl, file.path("results", "tables", "descriptive_statistics.csv"), row.names = FALSE)

# Compact national monthly table for the report
nat_month <- desc_tbl[desc_tbl$series %in% paste("India", MONTH_COLS), ]
write.csv(nat_month, file.path("results", "tables", "india_monthly_descriptive.csv"), row.names = FALSE)

nat_season <- desc_tbl[desc_tbl$series %in% c(
  "India annual (1901-2015)",
  "India Jan-Feb winter",
  "India Mar-May pre-monsoon",
  "India Jun-Sep southwest monsoon",
  "India Oct-Dec post-monsoon"
), ]
write.csv(nat_season, file.path("results", "tables", "india_seasonal_annual_descriptive.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Empirical probabilities (historical relative frequencies)
# ---------------------------------------------------------------------------
ann <- india$ANNUAL
jjas <- india$`Jun-Sep`
jf <- india$`Jan-Feb`
ond <- india$`Oct-Dec`
n <- length(ann)
lpa_ann <- mean(ann)
lpa_jjas <- mean(jjas)
p10 <- as.numeric(stats::quantile(ann, 0.10, names = FALSE))
p90 <- as.numeric(stats::quantile(ann, 0.90, names = FALSE))
p25_jf <- as.numeric(stats::quantile(jf, 0.25, names = FALSE))

cat_jjas <- monsoon_category(jjas, lpa_jjas)
cat_ann <- monsoon_category(ann, lpa_ann)

dry_winter <- jf < p25_jf
def_monsoon <- cat_jjas == "deficient"

prob_tbl <- rbind(
  data.frame(
    event = "Annual rainfall greater than the 1901-2015 sample mean (LPA analogue)",
    definition = paste0("ANNUAL > ", fmt(lpa_ann, 2), " mm"),
    numerator = sum(ann > lpa_ann),
    denominator = n,
    empirical_probability = mean(ann > lpa_ann),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Unusually dry year (below historical 10th percentile)",
    definition = paste0("ANNUAL < ", fmt(p10, 2), " mm (10th percentile of 1901-2015)"),
    numerator = sum(ann < p10),
    denominator = n,
    empirical_probability = mean(ann < p10),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Unusually wet year (above historical 90th percentile)",
    definition = paste0("ANNUAL > ", fmt(p90, 2), " mm (90th percentile of 1901-2015)"),
    numerator = sum(ann > p90),
    denominator = n,
    empirical_probability = mean(ann > p90),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Deficient southwest monsoon (IMD-style: JJAS < 90% of series mean)",
    definition = paste0("Jun-Sep < 0.90 * ", fmt(lpa_jjas, 2), " mm"),
    numerator = sum(cat_jjas == "deficient"),
    denominator = n,
    empirical_probability = mean(cat_jjas == "deficient"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Normal southwest monsoon (90-110% of series mean)",
    definition = "0.90*LPA <= Jun-Sep <= 1.10*LPA",
    numerator = sum(cat_jjas == "normal"),
    denominator = n,
    empirical_probability = mean(cat_jjas == "normal"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Excess southwest monsoon (JJAS > 110% of series mean)",
    definition = paste0("Jun-Sep > 1.10 * ", fmt(lpa_jjas, 2), " mm"),
    numerator = sum(cat_jjas == "excess"),
    denominator = n,
    empirical_probability = mean(cat_jjas == "excess"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Conditional: deficient monsoon given dry winter (JF < 25th percentile)",
    definition = paste0("P(JJAS deficient | Jan-Feb < ", fmt(p25_jf, 2), " mm)"),
    numerator = sum(def_monsoon & dry_winter),
    denominator = sum(dry_winter),
    empirical_probability = mean(def_monsoon[dry_winter]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Conditional: deficient monsoon given not-dry winter",
    definition = "P(JJAS deficient | Jan-Feb >= 25th percentile)",
    numerator = sum(def_monsoon & !dry_winter),
    denominator = sum(!dry_winter),
    empirical_probability = mean(def_monsoon[!dry_winter]),
    stringsAsFactors = FALSE
  ),
  data.frame(
    event = "Southwest monsoon is the wettest IMD season of the year",
    definition = "Jun-Sep == max(Jan-Feb, Mar-May, Jun-Sep, Oct-Dec)",
    numerator = sum(jjas >= india$`Jan-Feb` & jjas >= india$`Mar-May` & jjas >= ond),
    denominator = n,
    empirical_probability = mean(jjas >= india$`Jan-Feb` & jjas >= india$`Mar-May` & jjas >= ond),
    stringsAsFactors = FALSE
  )
)
write.csv(prob_tbl, file.path("results", "tables", "empirical_probabilities.csv"), row.names = FALSE)

freq_monsoon <- as.data.frame(table(factor(cat_jjas, levels = c("deficient", "normal", "excess"))),
                              stringsAsFactors = FALSE)
names(freq_monsoon) <- c("category", "count")
freq_monsoon$proportion <- freq_monsoon$count / n
freq_monsoon$lpa_mm <- lpa_jjas
write.csv(freq_monsoon, file.path("results", "tables", "monsoon_category_frequencies.csv"), row.names = FALSE)

# Frequency of dry months in arid vs wet subdivisions (zeros are genuine)
zero_focus <- c("West Rajasthan", "Kerala", "Tamil Nadu", "Konkan & Goa")
zero_rows <- list()
for (rg in zero_focus) {
  d <- sub_ok[sub_ok$SUBDIVISION == rg, ]
  for (m in MONTH_COLS) {
    x <- d[[m]]
    nobs <- sum(is.finite(x))
    n0 <- sum(is.finite(x) & x == 0)
    zero_rows[[length(zero_rows) + 1L]] <- data.frame(
      SUBDIVISION = rg, month = m, n_complete_years = nobs,
      n_zero = n0, empirical_p_zero = ifelse(nobs, n0 / nobs, NA_real_),
      stringsAsFactors = FALSE
    )
  }
}
zero_tbl <- do.call(rbind, zero_rows)
write.csv(zero_tbl, file.path("results", "tables", "empirical_p_zero_rainfall_selected_regions.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Commentary file: what the measures mean
# ---------------------------------------------------------------------------
ann_d <- desc_tbl[desc_tbl$series == "India annual (1901-2015)", ]
jjas_d <- desc_tbl[desc_tbl$series == "India Jun-Sep southwest monsoon", ]
write_txt(file.path("results", "statistical_outputs", "descriptive_and_probability_notes.txt"), c(
  "DESCRIPTIVE STATISTICS AND EMPIRICAL PROBABILITIES",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Units are millimetres. National series: 115 complete years, 1901-2015.",
  "Subdivision summaries use years with all 12 months observed, 1901-2015.",
  "",
  "INTERPRETATION OF MEASURES",
  "Mean: centre of mass of the historical series; sensitive to wet extremes.",
  "Median: 50% of years are wetter, 50% drier; more robust than the mean.",
  "Mode (rounded to 0.1 mm): often undefined for a continuous series with unique values;",
  "  NA in the table means no unique most-frequent rounded value.",
  "Range / min / max: observed envelope, not a prediction of future extremes.",
  "Variance and SD: year-to-year (or month-to-month) spread in mm (SD) or mm^2 (variance).",
  "Quartiles, IQR, percentiles: describe the bulk of the distribution without assuming normality.",
  "CV = 100*SD/mean: relative variability, comparable across wet and dry regions.",
  "Skewness: rainfall is typically right-skewed (long wet tail).",
  "Excess kurtosis: tail weight relative to a normal distribution (normal = 0).",
  "SE = SD/sqrt(n): sampling variability of the mean if years were iid; see script 05",
  "  for why that assumption is only approximate in a time series.",
  "",
  paste0("India annual mean = ", fmt(ann_d$mean, 2), " mm; median = ", fmt(ann_d$median, 2),
         " mm; SD = ", fmt(ann_d$sd, 2), " mm; CV = ", fmt(ann_d$cv_pct, 2), "%; skewness = ",
         fmt(ann_d$skewness, 3), "."),
  paste0("India Jun-Sep mean = ", fmt(jjas_d$mean, 2), " mm (about ",
         fmt(100 * jjas_d$mean / ann_d$mean, 1), "% of mean annual), SD = ",
         fmt(jjas_d$sd, 2), " mm."),
  "",
  "EMPIRICAL PROBABILITIES",
  "These are historical relative frequencies in 1901-2015, not proven future probabilities.",
  "They would equal long-run probabilities only under stationarity and iid sampling,",
  "which rainfall time series need not satisfy (see time-series script).",
  "IMD-style deficient/normal/excess uses 90% and 110% of THIS series mean as LPA analogue,",
  "not the official IMD 50-year long-period average, which is not stored in these files.",
  "",
  paste0("P(deficient JJAS) = ", fmt(mean(cat_jjas == "deficient"), 3),
         "; P(normal) = ", fmt(mean(cat_jjas == "normal"), 3),
         "; P(excess) = ", fmt(mean(cat_jjas == "excess"), 3), "."),
  paste0("P(deficient monsoon | dry winter) = ",
         fmt(mean(def_monsoon[dry_winter]), 3),
         " versus P(deficient | not dry winter) = ",
         fmt(mean(def_monsoon[!dry_winter]), 3), "."),
  "A difference, if any, is descriptive association, not a causal claim about winter rain",
  "controlling the monsoon."
))

# Monthly mean barplot
save_png("india_monthly_mean_rainfall.png", {
  means <- sapply(MONTH_COLS, function(m) mean(india[[m]]))
  graphics::barplot(means, names.arg = MONTH_LABELS, col = "#3B6D9A", border = NA,
                    ylab = "Mean rainfall (mm)", xlab = "Month",
                    main = "India area-weighted mean monthly rainfall, 1901-2015")
  graphics::abline(h = mean(means), lty = 2, col = "grey30")
})

save_png("india_seasonal_boxplot.png", {
  seas <- data.frame(
    season = factor(rep(c("Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec"), each = n),
                    levels = c("Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")),
    rainfall = c(india$`Jan-Feb`, india$`Mar-May`, india$`Jun-Sep`, india$`Oct-Dec`)
  )
  graphics::boxplot(rainfall ~ season, data = seas, col = c("#A6CEE3", "#B2DF8A", "#1F78B4", "#FDBF6F"),
                    ylab = "Rainfall (mm)", xlab = "IMD season",
                    main = "India seasonal rainfall distributions, 1901-2015")
})

save_png("india_monthly_boxplots.png", {
  mdat <- as.list(india[MONTH_COLS])
  names(mdat) <- MONTH_LABELS
  graphics::boxplot(mdat, col = "#9ECAE1", ylab = "Rainfall (mm)", xlab = "Month",
                    main = "India monthly rainfall boxplots, 1901-2015", las = 2)
})

cat("03_descriptive_analysis.R complete.\n")
