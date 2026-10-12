# =============================================================================
# 02_data_cleaning.R
# Phase 1 cleaning: analysis copies only. Original CSVs are never overwritten.
#
# Prerequisites: 01_data_inspection.R
# Run from project root:
#   Rscript R_code/02_data_cleaning.R
# =============================================================================

set_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    setwd(dirname(dirname(script_path)))
    return(invisible(getwd()))
  }
  if (basename(getwd()) == "R_code") setwd("..")
  invisible(getwd())
}

set_project_root()

na_strings <- c("NA", "NaN", "N/A", "n/a", "NULL", "null", "-", "--", ".", "")
month_cols <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

season_map_national <- list(
  `Jan-Feb` = c("JAN", "FEB"),
  `Mar-May` = c("MAR", "APR", "MAY"),
  `Jun-Sep` = c("JUN", "JUL", "AUG", "SEP"),
  `Oct-Dec` = c("OCT", "NOV", "DEC")
)
season_map_subdiv <- list(
  JF   = c("JAN", "FEB"),
  MAM  = c("MAR", "APR", "MAY"),
  JJAS = c("JUN", "JUL", "AUG", "SEP"),
  OND  = c("OCT", "NOV", "DEC")
)

ensure_dirs <- function() {
  dirs <- c(file.path("results", "tables"),
            file.path("results", "statistical_outputs"),
            file.path("results", "cleaned_data"))
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
ensure_dirs()

read_rainfall_csv <- function(path) {
  read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = na_strings,
    strip.white = TRUE
  )
}

cleaning_log <- data.frame(
  step_id = integer(),
  dataset = character(),
  action = character(),
  justification = character(),
  n_affected = integer(),
  details = character(),
  stringsAsFactors = FALSE
)

append_log <- function(step_id, dataset, action, justification, n_affected, details = "") {
  cleaning_log <<- rbind(
    cleaning_log,
    data.frame(
      step_id = as.integer(step_id),
      dataset = dataset,
      action = action,
      justification = justification,
      n_affected = as.integer(n_affected),
      details = details,
      stringsAsFactors = FALSE
    )
  )
}

to_numeric_rainfall <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

row_sums_checked <- function(df, cols) {
  m <- as.data.frame(lapply(df[cols], to_numeric_rainfall), stringsAsFactors = FALSE)
  rowSums(m, na.rm = FALSE)
}

iqr_outlier_flag <- function(x, k = 3) {
  xnum <- to_numeric_rainfall(x)
  out <- rep(FALSE, length(xnum))
  ok <- is.finite(xnum)
  if (sum(ok) < 12L) return(out)
  qs <- stats::quantile(xnum[ok], probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE, type = 7)
  iqr <- qs[2] - qs[1]
  if (!is.finite(iqr) || iqr == 0) {
    # Fall back to a robust scale so constant series are not all flagged.
    s <- stats::mad(xnum[ok], constant = 1)
    if (!is.finite(s) || s == 0) return(out)
    med <- stats::median(xnum[ok])
    out[ok] <- abs(xnum[ok] - med) > (k * 6 * s / 1.349)
    return(out)
  }
  lo <- qs[1] - k * iqr
  hi <- qs[2] + k * iqr
  out[ok] <- xnum[ok] < lo | xnum[ok] > hi
  out
}

# ---------------------------------------------------------------------------
# National series
# ---------------------------------------------------------------------------

national_path <- file.path("data", "rainfall_area-wt_India_1901-2015.csv")
india_raw <- read_rainfall_csv(national_path)
n0 <- nrow(india_raw)
append_log(1, "national_area_weighted", "import",
           "Read original CSV with explicit NA strings; originals not modified.",
           n0, paste("file=", national_path))

names(india_raw) <- trimws(names(india_raw))
india <- india_raw

india$REGION <- trimws(as.character(india$REGION))
n_ws_region <- sum(india_raw$REGION != india$REGION, na.rm = TRUE)
append_log(2, "national_area_weighted", "trim_REGION",
           "Remove leading/trailing whitespace from region labels.",
           n_ws_region, paste(unique(india$REGION), collapse = ", "))

year_before <- india$YEAR
india$YEAR <- as.integer(round(to_numeric_rainfall(india$YEAR)))
append_log(3, "national_area_weighted", "YEAR_to_integer",
           "Calendar year is a discrete index, not a continuous measurement.",
           sum(is.na(india$YEAR) != is.na(year_before) | (!is.na(india$YEAR) & !is.na(as.numeric(year_before)) & india$YEAR != as.integer(as.numeric(year_before)))),
           paste0("range=", min(india$YEAR, na.rm = TRUE), "-", max(india$YEAR, na.rm = TRUE)))

rain_cols_nat <- c(month_cols, "ANNUAL", "Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")
n_coerced <- 0L
for (col in rain_cols_nat) {
  if (!col %in% names(india)) next
  before <- india[[col]]
  india[[col]] <- to_numeric_rainfall(before)
  n_coerced <- n_coerced + sum(is.na(india[[col]]) & !is.na(before))
}
append_log(4, "national_area_weighted", "rainfall_to_numeric",
           "Rainfall is a continuous millimetre measurement; coercion failures would become NA (none expected).",
           n_coerced, paste(rain_cols_nat, collapse = ", "))

n_neg <- sum(sapply(india[intersect(month_cols, names(india))], function(z) sum(!is.na(z) & z < 0)))
append_log(5, "national_area_weighted", "negative_rainfall_screen",
           "Negative millimetres are physically impossible; none are deleted because none were found, or they would be flagged not dropped.",
           n_neg, "threshold: month < 0 mm")

dup_exact <- duplicated(india) | duplicated(india, fromLast = TRUE)
dup_key <- duplicated(india[, c("REGION", "YEAR")]) | duplicated(india[, c("REGION", "YEAR")], fromLast = TRUE)
append_log(6, "national_area_weighted", "duplicate_investigation",
           "Exact row duplicates and REGION-YEAR duplicates were counted. No rows removed because none were found.",
           sum(dup_key), paste0("exact_dup_rows=", sum(duplicated(india)), "; key_dups=", sum(duplicated(india[, c("REGION", "YEAR")]))))

india$n_missing_months <- rowSums(is.na(india[month_cols]))
india$annual_computed <- row_sums_checked(india, month_cols)
india$janfeb_computed <- row_sums_checked(india, c("JAN", "FEB"))
india$marmay_computed <- row_sums_checked(india, c("MAR", "APR", "MAY"))
india$junsep_computed <- row_sums_checked(india, c("JUN", "JUL", "AUG", "SEP"))
india$octdec_computed <- row_sums_checked(india, c("OCT", "NOV", "DEC"))

india$flag_inconsistent_annual <- with(india, {
  ok <- is.finite(ANNUAL) & is.finite(annual_computed)
  ifelse(ok, abs(ANNUAL - annual_computed) > 1, FALSE)
})
india$flag_inconsistent_seasonal <- with(india, {
  d1 <- is.finite(`Jan-Feb`) & is.finite(janfeb_computed) & abs(`Jan-Feb` - janfeb_computed) > 1
  d2 <- is.finite(`Mar-May`) & is.finite(marmay_computed) & abs(`Mar-May` - marmay_computed) > 1
  d3 <- is.finite(`Jun-Sep`) & is.finite(junsep_computed) & abs(`Jun-Sep` - junsep_computed) > 1
  d4 <- is.finite(`Oct-Dec`) & is.finite(octdec_computed) & abs(`Oct-Dec` - octdec_computed) > 1
  d1 | d2 | d3 | d4
})

n_ann_inc <- sum(india$flag_inconsistent_annual, na.rm = TRUE)
n_sea_inc <- sum(india$flag_inconsistent_seasonal, na.rm = TRUE)
append_log(7, "national_area_weighted", "consistency_flags",
           "Source seasonal/annual totals were retained. Differences > 1 mm versus the sum of constituent months are flagged, not overwritten.",
           n_ann_inc + n_sea_inc,
           paste0("annual_inconsistent=", n_ann_inc, "; seasonal_inconsistent=", n_sea_inc, "; tolerance=1 mm"))

# Reconstruction only where source derived field is NA and months are complete.
india$annual_filled_from_months <- is.na(india$ANNUAL) & india$n_missing_months == 0L
india$ANNUAL_analysis <- ifelse(is.na(india$ANNUAL) & india$n_missing_months == 0L,
                                india$annual_computed, india$ANNUAL)
append_log(8, "national_area_weighted", "annual_reconstruct_if_source_NA",
           "Not imputation of unobserved rain: ANNUAL is a total of observed months. Filled only when all 12 months exist and source ANNUAL is NA.",
           sum(india$annual_filled_from_months),
           "New column ANNUAL_analysis; source ANNUAL unchanged.")

# Tukey fences by month across the national series (k=3, conservative).
india$n_month_outlier_flags <- 0L
for (m in month_cols) {
  fl <- iqr_outlier_flag(india[[m]], k = 3)
  india$n_month_outlier_flags <- india$n_month_outlier_flags + as.integer(fl)
}
append_log(9, "national_area_weighted", "outlier_flag_iqr_k3",
           "Extreme rainfall years are scientifically informative (floods/droughts). Flagged with a conservative 3xIQR rule; no rows deleted.",
           sum(india$n_month_outlier_flags > 0),
           paste0("rows_with_any_month_flag=", sum(india$n_month_outlier_flags > 0),
                  "; total_month_flags=", sum(india$n_month_outlier_flags)))

india$qc_status <- ifelse(
  india$n_missing_months == 0L & !india$flag_inconsistent_annual & !india$flag_inconsistent_seasonal,
  "complete_consistent",
  ifelse(india$n_missing_months > 0L, "incomplete_months", "inconsistent_totals")
)

india_out <- india[order(india$YEAR), ]
rownames(india_out) <- NULL
nat_clean_path <- file.path("results", "cleaned_data", "india_national_rainfall_cleaned.csv")
write.csv(india_out, nat_clean_path, row.names = FALSE)
append_log(10, "national_area_weighted", "write_cleaned_copy",
           "Analysis copy written; original data/ CSV untouched.",
           nrow(india_out), nat_clean_path)

# ---------------------------------------------------------------------------
# Subdivision series
# ---------------------------------------------------------------------------

subdiv_path <- file.path("data", "Sub_Division_IMD_2017.csv")
sub_raw <- read_rainfall_csv(subdiv_path)
append_log(11, "imd_subdivision", "import",
           "Read original CSV with explicit NA strings; originals not modified.",
           nrow(sub_raw), paste("file=", subdiv_path))

names(sub_raw) <- trimws(names(sub_raw))
subdiv <- sub_raw

label_before <- as.character(subdiv$SUBDIVISION)
subdiv$SUBDIVISION <- trimws(gsub("[[:space:]]+", " ", as.character(subdiv$SUBDIVISION)))
n_label <- sum(label_before != subdiv$SUBDIVISION, na.rm = TRUE)
append_log(12, "imd_subdivision", "standardise_SUBDIVISION_whitespace",
           "Trim and collapse repeated whitespace. Historical IMD names (Orissa, Matathwada, Rayalseema) are retained so records stay matchable to the source.",
           n_label, paste0("n_unique_labels=", length(unique(subdiv$SUBDIVISION))))

subdiv$YEAR <- as.integer(round(to_numeric_rainfall(subdiv$YEAR)))
append_log(13, "imd_subdivision", "YEAR_to_integer",
           "Calendar year is a discrete index.",
           sum(is.na(subdiv$YEAR)),
           paste0("range=", min(subdiv$YEAR, na.rm = TRUE), "-", max(subdiv$YEAR, na.rm = TRUE)))

rain_cols_sub <- c(month_cols, "ANNUAL", "JF", "MAM", "JJAS", "OND")
n_coerced_s <- 0L
for (col in rain_cols_sub) {
  before <- subdiv[[col]]
  subdiv[[col]] <- to_numeric_rainfall(before)
  n_coerced_s <- n_coerced_s + sum(is.na(subdiv[[col]]) & !is.na(before))
}
append_log(14, "imd_subdivision", "rainfall_to_numeric",
           "Rainfall coerced to numeric millimetres. Source 'NA' tokens already mapped to R NA on import.",
           n_coerced_s, paste(rain_cols_sub, collapse = ", "))

n_zero <- sum(sapply(subdiv[month_cols], function(z) sum(!is.na(z) & z == 0)))
n_miss_m <- sum(sapply(subdiv[month_cols], function(z) sum(is.na(z))))
append_log(15, "imd_subdivision", "retain_zeros_and_missing",
           "Zero millimetres remain as dry months. NA remains missing. Incomplete rows were not deleted.",
           n_miss_m, paste0("monthly_zeros=", n_zero, "; monthly_NAs=", n_miss_m,
                            "; rows_kept=", nrow(subdiv)))

n_neg_s <- sum(sapply(subdiv[month_cols], function(z) sum(!is.na(z) & z < 0)))
append_log(16, "imd_subdivision", "negative_rainfall_screen",
           "Negative rainfall would be flagged, not deleted.",
           n_neg_s, "threshold: month < 0 mm")

dup_exact_s <- sum(duplicated(subdiv[, c("SUBDIVISION", "YEAR", month_cols)]))
key <- subdiv[, c("SUBDIVISION", "YEAR")]
dup_key_s <- duplicated(key) | duplicated(key, fromLast = TRUE)
n_dup_keys <- sum(duplicated(key))
if (n_dup_keys > 0) {
  # Keep first occurrence only if entire rows are identical; otherwise keep all and flag.
  identical_dups <- duplicated(subdiv)
  if (all(duplicated(subdiv) | !dup_key_s[seq_len(nrow(subdiv))])) {
    n_before <- nrow(subdiv)
    subdiv <- subdiv[!duplicated(subdiv), ]
    append_log(17, "imd_subdivision", "drop_identical_duplicate_rows",
               "Removed extra copies of fully identical rows after confirming they were exact duplicates.",
               n_before - nrow(subdiv), "kept first occurrence")
  } else {
    subdiv$flag_duplicate_key <- dup_key_s
    append_log(17, "imd_subdivision", "flag_conflicting_duplicate_keys",
               "SUBDIVISION-YEAR repeats were not deleted because values may conflict; flagged for review.",
               sum(dup_key_s), paste0("n_duplicate_keys=", n_dup_keys))
  }
} else {
  subdiv$flag_duplicate_key <- FALSE
  append_log(17, "imd_subdivision", "duplicate_investigation",
             "No SUBDIVISION-YEAR duplicates. No rows removed.",
             0L, paste0("exact_key_dups=", n_dup_keys))
}

subdiv$n_missing_months <- rowSums(is.na(subdiv[month_cols]))
subdiv$annual_computed <- row_sums_checked(subdiv, month_cols)
subdiv$jf_computed   <- row_sums_checked(subdiv, season_map_subdiv$JF)
subdiv$mam_computed  <- row_sums_checked(subdiv, season_map_subdiv$MAM)
subdiv$jjas_computed <- row_sums_checked(subdiv, season_map_subdiv$JJAS)
subdiv$ond_computed  <- row_sums_checked(subdiv, season_map_subdiv$OND)

tol <- 1
subdiv$flag_inconsistent_annual <- with(subdiv, {
  ok <- is.finite(ANNUAL) & is.finite(annual_computed)
  ifelse(ok, abs(ANNUAL - annual_computed) > tol, FALSE)
})
subdiv$flag_inconsistent_seasonal <- with(subdiv, {
  d1 <- is.finite(JF) & is.finite(jf_computed) & abs(JF - jf_computed) > tol
  d2 <- is.finite(MAM) & is.finite(mam_computed) & abs(MAM - mam_computed) > tol
  d3 <- is.finite(JJAS) & is.finite(jjas_computed) & abs(JJAS - jjas_computed) > tol
  d4 <- is.finite(OND) & is.finite(ond_computed) & abs(OND - ond_computed) > tol
  d1 | d2 | d3 | d4
})
append_log(18, "imd_subdivision", "consistency_flags",
           "Source derived totals retained. Inconsistencies > 1 mm versus complete month sums are flagged, not overwritten.",
           sum(subdiv$flag_inconsistent_annual | subdiv$flag_inconsistent_seasonal, na.rm = TRUE),
           paste0("annual=", sum(subdiv$flag_inconsistent_annual, na.rm = TRUE),
                  "; seasonal=", sum(subdiv$flag_inconsistent_seasonal, na.rm = TRUE)))

# Fill analysis derived fields from months only when source is NA and parts are complete.
fill_from_months <- function(src, computed, needed_cols, df) {
  complete_parts <- rowSums(is.na(df[needed_cols])) == 0L
  ifelse(is.na(src) & complete_parts, computed, src)
}

subdiv$ANNUAL_analysis <- fill_from_months(subdiv$ANNUAL, subdiv$annual_computed, month_cols, subdiv)
subdiv$JF_analysis     <- fill_from_months(subdiv$JF, subdiv$jf_computed, season_map_subdiv$JF, subdiv)
subdiv$MAM_analysis    <- fill_from_months(subdiv$MAM, subdiv$mam_computed, season_map_subdiv$MAM, subdiv)
subdiv$JJAS_analysis   <- fill_from_months(subdiv$JJAS, subdiv$jjas_computed, season_map_subdiv$JJAS, subdiv)
subdiv$OND_analysis    <- fill_from_months(subdiv$OND, subdiv$ond_computed, season_map_subdiv$OND, subdiv)

subdiv$annual_filled_from_months <- is.na(subdiv$ANNUAL) & subdiv$n_missing_months == 0L
subdiv$jf_filled_from_months   <- is.na(subdiv$JF)   & rowSums(is.na(subdiv[season_map_subdiv$JF])) == 0L
subdiv$mam_filled_from_months  <- is.na(subdiv$MAM)  & rowSums(is.na(subdiv[season_map_subdiv$MAM])) == 0L
subdiv$jjas_filled_from_months <- is.na(subdiv$JJAS) & rowSums(is.na(subdiv[season_map_subdiv$JJAS])) == 0L
subdiv$ond_filled_from_months  <- is.na(subdiv$OND)  & rowSums(is.na(subdiv[season_map_subdiv$OND])) == 0L

append_log(19, "imd_subdivision", "reconstruct_derived_totals_if_source_NA",
           "Where IMD omitted ANNUAL/JF/MAM/JJAS/OND but all constituent months are observed, the analysis column is the month sum. This does not invent monthly rain.",
           sum(subdiv$annual_filled_from_months | subdiv$jf_filled_from_months |
                 subdiv$mam_filled_from_months | subdiv$jjas_filled_from_months |
                 subdiv$ond_filled_from_months),
           paste0("ANNUAL_filled=", sum(subdiv$annual_filled_from_months),
                  "; JF=", sum(subdiv$jf_filled_from_months),
                  "; MAM=", sum(subdiv$mam_filled_from_months),
                  "; JJAS=", sum(subdiv$jjas_filled_from_months),
                  "; OND=", sum(subdiv$ond_filled_from_months)))

# Outliers within each subdivision-month (k=3 IQR). Keep all rows.
subdiv$n_month_outlier_flags <- 0L
regs <- unique(subdiv$SUBDIVISION)
for (m in month_cols) {
  flags_m <- rep(FALSE, nrow(subdiv))
  for (rg in regs) {
    idx <- which(subdiv$SUBDIVISION == rg)
    flags_m[idx] <- iqr_outlier_flag(subdiv[[m]][idx], k = 3)
  }
  subdiv$n_month_outlier_flags <- subdiv$n_month_outlier_flags + as.integer(flags_m)
}
append_log(20, "imd_subdivision", "outlier_flag_iqr_k3_by_subdivision_month",
           "Outliers flagged within each subdivision and month using 3xIQR. Extreme monsoon totals are retained.",
           sum(subdiv$n_month_outlier_flags > 0),
           paste0("rows_with_any_month_flag=", sum(subdiv$n_month_outlier_flags > 0),
                  "; total_month_flags=", sum(subdiv$n_month_outlier_flags)))

# Unbalanced panel: document missing years but do not fabricate rows.
year_gaps <- 0L
gap_examples <- character(0)
for (rg in regs) {
  y <- sort(unique(subdiv$YEAR[subdiv$SUBDIVISION == rg]))
  y <- y[!is.na(y)]
  if (length(y) < 2L) next
  missing <- setdiff(seq(min(y), max(y)), y)
  year_gaps <- year_gaps + length(missing)
  if (length(missing) && length(gap_examples) < 8) {
    gap_examples <- c(gap_examples, paste0(rg, ": ", paste(missing, collapse = ",")))
  }
}
append_log(21, "imd_subdivision", "retain_unbalanced_panel",
           "Internal missing years (e.g. a subdivision with 1908 then 1910) are left as absent rows. No years were invented.",
           year_gaps, paste(gap_examples, collapse = " | "))

subdiv$months_complete <- subdiv$n_missing_months == 0L
subdiv$qc_status <- ifelse(
  subdiv$n_missing_months == 12L,
  "all_months_missing",
  ifelse(
    subdiv$n_missing_months > 0L,
    "incomplete_months",
    ifelse(subdiv$flag_inconsistent_annual | subdiv$flag_inconsistent_seasonal,
           "complete_months_inconsistent_totals",
           "complete_consistent")
  )
)

sub_out <- subdiv[order(subdiv$SUBDIVISION, subdiv$YEAR), ]
rownames(sub_out) <- NULL
sub_clean_path <- file.path("results", "cleaned_data", "imd_subdivision_rainfall_cleaned.csv")
write.csv(sub_out, sub_clean_path, row.names = FALSE)
append_log(22, "imd_subdivision", "write_cleaned_copy",
           "Analysis copy written; original data/ CSV untouched.",
           nrow(sub_out), sub_clean_path)

# ---------------------------------------------------------------------------
# Persist log and append to quality summary
# ---------------------------------------------------------------------------

log_path <- file.path("results", "tables", "cleaning_log.csv")
write.csv(cleaning_log, log_path, row.names = FALSE)

qc_path <- file.path("results", "statistical_outputs", "data_quality_summary.txt")
sink(qc_path, append = TRUE)
cat("\n==============================================================================\n")
cat("6. CLEANING DECISIONS AND COUNTS (from 02_data_cleaning.R)\n")
cat("==============================================================================\n")
cat("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
print(cleaning_log, right = FALSE, row.names = FALSE)
cat("\nCleaned analysis files:\n")
cat("  ", nat_clean_path, "  rows=", nrow(india_out), " cols=", ncol(india_out), "\n", sep = "")
cat("  ", sub_clean_path, "  rows=", nrow(sub_out), " cols=", ncol(sub_out), "\n", sep = "")
cat("\nUnresolved issues carried into later phases:\n")
cat("  - Subdivision panel is unbalanced (late starts and internal year gaps).\n")
cat("  - Monthly NA values remain NA; analysis of annual/seasonal totals should\n")
cat("    use *_analysis columns and/or restrict to months_complete == TRUE.\n")
cat("  - A small number of source seasonal/annual totals may differ from month\n")
cat("    sums by more than 1 mm; source values are kept and flagged.\n")
cat("  - National and subdivision series have different spatial support and\n")
cat("    slightly different end years (2015 vs 2017).\n")
cat("  - Original files in data/ were not modified.\n")
sink()

cat("Cleaning complete.\n")
cat("Wrote: ", log_path, "\n", sep = "")
cat("Wrote: ", nat_clean_path, "\n", sep = "")
cat("Wrote: ", sub_clean_path, "\n", sep = "")
cat("Updated: ", qc_path, "\n", sep = "")
