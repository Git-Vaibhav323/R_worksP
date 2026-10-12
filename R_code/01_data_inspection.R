# =============================================================================
# 01_data_inspection.R
# Statistical Analysis of Rainfall Variability and Regional Monsoon Patterns
# in India — Phase 1 data-quality audit (no modelling).
#
# Run from the project root in RStudio (Session > Set Working Directory >
# To Source File Location, then setwd("..") if needed), or:
#   Rscript R_code/01_data_inspection.R
# =============================================================================

set_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    setwd(dirname(dirname(script_path)))
    return(invisible(getwd()))
  }
  if (basename(getwd()) == "R_code") {
    setwd("..")
  }
  invisible(getwd())
}

set_project_root()

na_strings <- c("NA", "NaN", "N/A", "n/a", "NULL", "null", "-", "--", ".", "")

month_cols <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN",
                "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")

ensure_dirs <- function() {
  dirs <- c("results",
            file.path("results", "tables"),
            file.path("results", "statistical_outputs"),
            file.path("results", "cleaned_data"),
            "graphs",
            "report",
            "R_code",
            "data")
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

ensure_dirs()

read_rainfall_csv <- function(path) {
  df <- read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = na_strings,
    strip.white = TRUE
  )
  names(df) <- trimws(names(df))
  df
}

is_numeric_like <- function(x) {
  if (is.numeric(x)) return(TRUE)
  xc <- trimws(as.character(x))
  xc <- xc[!is.na(xc) & nzchar(xc)]
  if (!length(xc)) return(FALSE)
  suppressWarnings(all(!is.na(as.numeric(xc))))
}

fmt_pct <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), 100 * x)
}

safe_range <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  range(x)
}

count_blank_raw <- function(path, col_index) {
  raw <- read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = NULL,
    colClasses = "character"
  )
  names(raw) <- trimws(names(raw))
  x <- raw[[col_index]]
  if (is.null(x)) return(0L)
  sum(is.na(x) | !nzchar(trimws(x)))
}

season_defs <- list(
  national = list(
    `Jan-Feb` = c("JAN", "FEB"),
    `Mar-May` = c("MAR", "APR", "MAY"),
    `Jun-Sep` = c("JUN", "JUL", "AUG", "SEP"),
    `Oct-Dec` = c("OCT", "NOV", "DEC"),
    ANNUAL   = month_cols
  ),
  subdivision = list(
    JF    = c("JAN", "FEB"),
    MAM   = c("MAR", "APR", "MAY"),
    JJAS  = c("JUN", "JUL", "AUG", "SEP"),
    OND   = c("OCT", "NOV", "DEC"),
    ANNUAL = month_cols
  )
)

check_derived_totals <- function(df, defs, tolerance_mm = 1) {
  out <- data.frame(
    derived_column = names(defs),
    n_comparable = NA_integer_,
    n_inconsistent = NA_integer_,
    max_abs_diff_mm = NA_real_,
    median_abs_diff_mm = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(defs)) {
    nm <- names(defs)[i]
    parts <- defs[[i]]
    if (!all(c(nm, parts) %in% names(df))) next
    src <- suppressWarnings(as.numeric(df[[nm]]))
    computed <- rowSums(as.data.frame(lapply(df[parts], function(z) {
      suppressWarnings(as.numeric(z))
    })), na.rm = FALSE)
    ok <- is.finite(src) & is.finite(computed)
    diffs <- abs(src[ok] - computed[ok])
    out$n_comparable[i] <- sum(ok)
    out$n_inconsistent[i] <- sum(diffs > tolerance_mm)
    out$max_abs_diff_mm[i] <- if (length(diffs)) max(diffs) else NA_real_
    out$median_abs_diff_mm[i] <- if (length(diffs)) stats::median(diffs) else NA_real_
  }
  out
}

national_path <- file.path("data", "rainfall_area-wt_India_1901-2015.csv")
subdiv_path   <- file.path("data", "Sub_Division_IMD_2017.csv")

stopifnot(file.exists(national_path), file.exists(subdiv_path))

india <- read_rainfall_csv(national_path)
subdiv <- read_rainfall_csv(subdiv_path)

# ---------------------------------------------------------------------------
# Structure summaries
# ---------------------------------------------------------------------------

structure_one <- function(df, dataset, path) {
  n <- nrow(df)
  p <- ncol(df)
  rows <- vector("list", p)
  for (j in seq_len(p)) {
    col <- names(df)[j]
    x <- df[[j]]
    xnum <- suppressWarnings(as.numeric(x))
    rng <- safe_range(x)
    n_miss <- sum(is.na(x) | (is.character(x) & !nzchar(trimws(x))))
    n_unique <- length(unique(x[!is.na(x)]))
    example <- NA_character_
    nonmiss <- x[!is.na(x)]
    if (length(nonmiss)) example <- as.character(nonmiss[1])
    rows[[j]] <- data.frame(
      dataset = dataset,
      source_file = basename(path),
      n_rows = n,
      n_cols = p,
      column_position = j,
      column_name = col,
      r_class = paste(class(x), collapse = "/"),
      treated_as_numeric = is_numeric_like(x),
      n_missing = n_miss,
      pct_missing = round(100 * n_miss / n, 4),
      n_unique_nonmissing = n_unique,
      min_numeric = rng[1],
      max_numeric = rng[2],
      example_value = example,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

structure_tbl <- rbind(
  structure_one(india, "national_area_weighted", national_path),
  structure_one(subdiv, "imd_subdivision", subdiv_path)
)
write.csv(
  structure_tbl,
  file.path("results", "tables", "data_structure_summary.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Missing values, zeros, negatives, special codes
# ---------------------------------------------------------------------------

special_code_counts <- function(path) {
  raw <- read.csv(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = NULL,
    colClasses = "character"
  )
  names(raw) <- trimws(names(raw))
  codes <- c("NA", "NaN", "N/A", "n/a", "NULL", "-", "--", ".", "NAA", "-999", "-99")
  out <- list()
  k <- 1L
  for (col in names(raw)) {
    x <- trimws(raw[[col]])
    blanks <- sum(!nzchar(x) | is.na(raw[[col]]))
    if (blanks > 0) {
      out[[k]] <- data.frame(column = col, code = "<blank>", n = blanks, stringsAsFactors = FALSE)
      k <- k + 1L
    }
    for (cd in codes) {
      n <- sum(tolower(x) == tolower(cd), na.rm = TRUE)
      if (n > 0) {
        out[[k]] <- data.frame(column = col, code = cd, n = n, stringsAsFactors = FALSE)
        k <- k + 1L
      }
    }
  }
  if (!length(out)) {
    return(data.frame(column = character(), code = character(), n = integer(), stringsAsFactors = FALSE))
  }
  do.call(rbind, out)
}

missing_one <- function(df, dataset) {
  n <- nrow(df)
  rows <- vector("list", ncol(df))
  for (j in seq_len(ncol(df))) {
    col <- names(df)[j]
    x <- df[[j]]
    xnum <- suppressWarnings(as.numeric(as.character(x)))
    n_miss <- sum(is.na(x))
    n_zero <- sum(!is.na(xnum) & xnum == 0)
    n_neg  <- sum(!is.na(xnum) & xnum < 0)
    n_nonmiss <- n - n_miss
    rows[[j]] <- data.frame(
      dataset = dataset,
      column = col,
      n_rows = n,
      n_missing = n_miss,
      pct_missing = round(100 * n_miss / n, 4),
      n_zero = n_zero,
      pct_zero_among_nonmissing = ifelse(n_nonmiss > 0, round(100 * n_zero / n_nonmiss, 4), NA_real_),
      n_negative = n_neg,
      n_nonmissing = n_nonmiss,
      note = ifelse(n_zero > 0,
                    "Zero values are retained as genuine zero rainfall, not missing.",
                    ""),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

missing_tbl <- rbind(
  missing_one(india, "national_area_weighted"),
  missing_one(subdiv, "imd_subdivision")
)
write.csv(
  missing_tbl,
  file.path("results", "tables", "missing_values_summary.csv"),
  row.names = FALSE
)

india_codes <- special_code_counts(national_path)
subdiv_codes <- special_code_counts(subdiv_path)
if (nrow(india_codes)) india_codes$dataset <- "national_area_weighted"
if (nrow(subdiv_codes)) subdiv_codes$dataset <- "imd_subdivision"
special_tbl <- rbind(india_codes, subdiv_codes)
if (nrow(special_tbl)) {
  write.csv(
    special_tbl,
    file.path("results", "tables", "special_missing_codes.csv"),
    row.names = FALSE
  )
}

# ---------------------------------------------------------------------------
# Duplicates, years, regions, consistency, implausible values
# ---------------------------------------------------------------------------

dup_rows_india <- sum(duplicated(india))
dup_keys_india <- sum(duplicated(india[, c("REGION", "YEAR"), drop = FALSE]))
dup_rows_sub <- sum(duplicated(subdiv))
dup_keys_sub <- sum(duplicated(subdiv[, c("SUBDIVISION", "YEAR"), drop = FALSE]))

india_years <- sort(unique(suppressWarnings(as.numeric(india$YEAR))))
sub_years <- sort(unique(suppressWarnings(as.numeric(subdiv$YEAR))))

india_expected_years <- seq(min(india_years, na.rm = TRUE), max(india_years, na.rm = TRUE))
india_missing_years <- setdiff(india_expected_years, india_years)

sub_regions <- sort(unique(trimws(as.character(subdiv$SUBDIVISION))))
year_by_region <- aggregate(
  YEAR ~ SUBDIVISION,
  data = transform(subdiv, SUBDIVISION = trimws(as.character(SUBDIVISION)),
                   YEAR = suppressWarnings(as.numeric(YEAR))),
  FUN = function(y) {
    y <- sort(unique(y[!is.na(y)]))
    paste(c(
      paste0("n=", length(y)),
      paste0("min=", min(y)),
      paste0("max=", max(y)),
      paste0("gaps=", length(setdiff(seq(min(y), max(y)), y)))
    ), collapse = "; ")
  }
)

india_cons <- check_derived_totals(india, season_defs$national, tolerance_mm = 1)
sub_cons   <- check_derived_totals(subdiv, season_defs$subdivision, tolerance_mm = 1)
india_cons$dataset <- "national_area_weighted"
sub_cons$dataset <- "imd_subdivision"
cons_tbl <- rbind(india_cons, sub_cons)
write.csv(
  cons_tbl,
  file.path("results", "tables", "seasonal_annual_consistency.csv"),
  row.names = FALSE
)

implausible_flags <- function(df, dataset, region_col) {
  out <- list()
  k <- 1L
  for (m in month_cols) {
    if (!m %in% names(df)) next
    x <- suppressWarnings(as.numeric(df[[m]]))
    n_neg <- sum(!is.na(x) & x < 0)
    n_huge <- sum(!is.na(x) & x > 4000)
    if (n_neg + n_huge > 0) {
      out[[k]] <- data.frame(
        dataset = dataset,
        column = m,
        n_negative = n_neg,
        n_gt_4000mm = n_huge,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  if (!length(out)) {
    return(data.frame(
      dataset = dataset, column = NA_character_,
      n_negative = 0L, n_gt_4000mm = 0L,
      stringsAsFactors = FALSE
    )[0, ])
  }
  do.call(rbind, out)
}

impl_india <- implausible_flags(india, "national_area_weighted", "REGION")
impl_sub <- implausible_flags(subdiv, "imd_subdivision", "SUBDIVISION")

# Completeness of monthly fields
n_missing_months_india <- rowSums(is.na(india[month_cols]))
n_missing_months_sub <- rowSums(is.na(subdiv[intersect(month_cols, names(subdiv))]))

region_completeness <- aggregate(
  n_missing_months_sub,
  by = list(SUBDIVISION = trimws(as.character(subdiv$SUBDIVISION))),
  FUN = function(z) {
    c(n_rows = length(z),
      n_complete_years = sum(z == 0),
      n_partial_years = sum(z > 0 & z < 12),
      n_all_months_missing = sum(z == 12),
      mean_missing_months = mean(z))
  }
)
region_completeness <- do.call(data.frame, region_completeness)
names(region_completeness) <- c(
  "SUBDIVISION", "n_rows", "n_complete_years",
  "n_partial_years", "n_all_months_missing", "mean_missing_months"
)
write.csv(
  region_completeness,
  file.path("results", "tables", "subdivision_year_completeness.csv"),
  row.names = FALSE
)

write.csv(
  year_by_region,
  file.path("results", "tables", "subdivision_year_coverage.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Text quality report
# ---------------------------------------------------------------------------

report_path <- file.path("results", "statistical_outputs", "data_quality_summary.txt")
con <- file(report_path, open = "wt")
on.exit(close(con), add = TRUE)

w <- function(...) {
  cat(..., "\n", sep = "", file = con)
}

w("DATA QUALITY SUMMARY")
w("Project: Statistical Analysis of Rainfall Variability and Regional Monsoon Patterns in India")
w("Phase: 1 — inspection and audit (no statistical modelling)")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
w("Working directory: ", getwd())
w("R version: ", paste(R.version$major, R.version$minor, sep = "."))
w("")
w("==============================================================================")
w("1. SOURCE FILES (UNMODIFIED ORIGINALS)")
w("==============================================================================")
w("National series : ", national_path,
  "  [size ", file.info(national_path)$size, " bytes]")
w("Subdivision series: ", subdiv_path,
  "  [size ", file.info(subdiv_path)$size, " bytes]")
w("")
w("Provenance (from file names and schema, not assumed extra metadata):")
w("  - Area-weighted all-India monthly/seasonal/annual rainfall, 1901-2015.")
w("  - India Meteorological Department (IMD) meteorological subdivision rainfall,")
w("    file stamp 2017, years present in the file: ",
  min(sub_years, na.rm = TRUE), "-", max(sub_years, na.rm = TRUE), ".")
w("  Units are millimetres of rainfall (standard IMD reporting). Zero is a valid")
w("  measurement (no rain), distinct from NA (not reported / not available).")
w("")
w("==============================================================================")
w("2. NATIONAL AREA-WEIGHTED DATASET")
w("==============================================================================")
w("Dimensions: ", nrow(india), " rows x ", ncol(india), " columns")
w("Column names: ", paste(names(india), collapse = " | "))
w("Identifying fields: REGION (constant geographical coverage), YEAR")
w("Monthly rainfall: JAN..DEC")
w("Seasonal rainfall: Jan-Feb (winter), Mar-May (pre-monsoon),")
w("                   Jun-Sep (southwest monsoon / JJAS), Oct-Dec (post-monsoon)")
w("Annual rainfall: ANNUAL")
w("")
w("REGION values: ", paste(unique(trimws(india$REGION)), collapse = ", "))
w("YEAR range: ", min(india_years), " to ", max(india_years),
  "  (n unique years = ", length(india_years), ")")
w("Missing years inside the min-max span: ",
  if (length(india_missing_years)) paste(india_missing_years, collapse = ", ") else "none")
w("Duplicate rows: ", dup_rows_india)
w("Duplicate REGION-YEAR keys: ", dup_keys_india)
w("Rows with any missing month: ", sum(n_missing_months_india > 0))
w("Rows with all 12 months present: ", sum(n_missing_months_india == 0))
w("Missing values by column:")
for (j in seq_len(ncol(india))) {
  nm <- names(india)[j]
  nmiss <- sum(is.na(india[[j]]))
  w("  ", sprintf("%-12s", nm), " missing = ", nmiss)
}
w("")
w("Numeric ranges (non-missing):")
for (col in c(month_cols, "ANNUAL", "Jan-Feb", "Mar-May", "Jun-Sep", "Oct-Dec")) {
  if (!col %in% names(india)) next
  rng <- safe_range(india[[col]])
  n0 <- sum(!is.na(suppressWarnings(as.numeric(india[[col]]))) &
              as.numeric(india[[col]]) == 0)
  w("  ", sprintf("%-12s", col),
    " min = ", sprintf("%.2f", rng[1]),
    "  max = ", sprintf("%.2f", rng[2]),
    "  zeros = ", n0)
}
w("")
w("Seasonal/annual vs constituent months (tolerance 1.0 mm):")
for (i in seq_len(nrow(india_cons))) {
  w("  ", india_cons$derived_column[i],
    ": comparable = ", india_cons$n_comparable[i],
    "; inconsistent = ", india_cons$n_inconsistent[i],
    "; max |diff| = ", sprintf("%.4f", india_cons$max_abs_diff_mm[i]), " mm",
    "; median |diff| = ", sprintf("%.4f", india_cons$median_abs_diff_mm[i]), " mm")
}
w("")
w("==============================================================================")
w("3. IMD SUBDIVISION DATASET")
w("==============================================================================")
w("Dimensions: ", nrow(subdiv), " rows x ", ncol(subdiv), " columns")
w("Column names: ", paste(names(subdiv), collapse = " | "))
w("Identifying fields: SUBDIVISION, YEAR")
w("Monthly rainfall: JAN..DEC")
w("Seasonal rainfall: JF, MAM, JJAS, OND  (same IMD season grouping as national)")
w("Annual rainfall: ANNUAL")
w("")
w("Unique SUBDIVISION labels (after trimws): ", length(sub_regions))
for (rg in sub_regions) w("  - ", rg)
w("")
w("YEAR range: ", min(sub_years), " to ", max(sub_years),
  "  (n unique years = ", length(sub_years), ")")
w("Duplicate rows: ", dup_rows_sub)
w("Duplicate SUBDIVISION-YEAR keys: ", dup_keys_sub)
w("Rows with all 12 months present: ", sum(n_missing_months_sub == 0))
w("Rows with 1-11 months missing: ", sum(n_missing_months_sub > 0 & n_missing_months_sub < 12))
w("Rows with all 12 months missing: ", sum(n_missing_months_sub == 12))
w("")
w("Missing values by column:")
for (j in seq_len(ncol(subdiv))) {
  nm <- names(subdiv)[j]
  nmiss <- sum(is.na(subdiv[[j]]))
  pct <- 100 * nmiss / nrow(subdiv)
  w("  ", sprintf("%-12s", nm), " missing = ", nmiss,
    "  (", sprintf("%.2f", pct), "%)")
}
w("")
w("Genuine zeros vs missing (monthly fields):")
for (m in month_cols) {
  x <- suppressWarnings(as.numeric(subdiv[[m]]))
  w("  ", m, ": zeros = ", sum(!is.na(x) & x == 0),
    "; missing = ", sum(is.na(x)),
    "; negative = ", sum(!is.na(x) & x < 0))
}
w("")
w("Implausible-value screen (month < 0 mm or month > 4000 mm):")
if (!nrow(impl_sub) && !nrow(impl_india)) {
  w("  None detected in either dataset.")
} else {
  if (nrow(impl_india)) print_target <- impl_india else print_target <- impl_sub
  if (nrow(impl_india)) {
    for (i in seq_len(nrow(impl_india))) {
      w("  national ", impl_india$column[i],
        " negatives=", impl_india$n_negative[i],
        " >4000mm=", impl_india$n_gt_4000mm[i])
    }
  }
  if (nrow(impl_sub)) {
    for (i in seq_len(nrow(impl_sub))) {
      w("  subdivision ", impl_sub$column[i],
        " negatives=", impl_sub$n_negative[i],
        " >4000mm=", impl_sub$n_gt_4000mm[i])
    }
  }
}
w("")
w("Seasonal/annual vs constituent months (tolerance 1.0 mm; only rows where")
w("all constituent months AND the derived field are non-missing):")
for (i in seq_len(nrow(sub_cons))) {
  w("  ", sub_cons$derived_column[i],
    ": comparable = ", sub_cons$n_comparable[i],
    "; inconsistent = ", sub_cons$n_inconsistent[i],
    "; max |diff| = ", sprintf("%.4f", sub_cons$max_abs_diff_mm[i]), " mm",
    "; median |diff| = ", sprintf("%.4f", sub_cons$median_abs_diff_mm[i]), " mm")
}
w("")
w("Coverage is not a balanced panel. Several subdivisions start after 1901")
w("and/or contain internal year gaps (for example Andaman & Nicobar Islands")
w("has no 1909 row in the source file). See:")
w("  results/tables/subdivision_year_coverage.csv")
w("  results/tables/subdivision_year_completeness.csv")
w("")
w("==============================================================================")
w("4. CROSS-DATASET NOTES")
w("==============================================================================")
w("Geographical coverage:")
w("  National file REGION = INDIA (one all-India area-weighted series).")
w("  Subdivision file contains ", length(sub_regions), " IMD meteorological subdivisions.")
w("  These are NOT the same spatial support: subdivision values must not be")
w("  averaged with equal weights to reconstruct the national series.")
w("Time span:")
w("  National: ", min(india_years), "-", max(india_years), ".")
w("  Subdivision: ", min(sub_years), "-", max(sub_years), ".")
w("  Overlap for joint work: ", max(min(india_years), min(sub_years)), "-",
  min(max(india_years), max(sub_years)), ".")
w("Season labels differ in spelling (Jan-Feb vs JF) but match IMD definitions.")
w("Both files report rainfall in millimetres; monthly and seasonal fields are")
w("source-provided, not computed by this project, unless a later cleaning step")
w("explicitly reconstructs a derived total from complete months.")
w("")
w("==============================================================================")
w("5. CLEANING IMPLICATIONS (implemented in 02_data_cleaning.R)")
w("==============================================================================")
w("- Preserve original CSVs in data/.")
w("- Do not drop incomplete years; missing months remain NA.")
w("- Do not impute rainfall. Derived seasonal/annual fields may be reconstructed")
w("  from complete constituent months into NEW columns only.")
w("- Do not delete statistical outliers; flag them.")
w("- Trim region labels; keep historical IMD names (e.g. Orissa, Matathwada).")
w("- Treat zeros as valid dry-month observations.")
w("")
w("Inspection artefacts:")
w("  results/tables/data_structure_summary.csv")
w("  results/tables/missing_values_summary.csv")
w("  results/tables/seasonal_annual_consistency.csv")
w("  results/tables/subdivision_year_coverage.csv")
w("  results/tables/subdivision_year_completeness.csv")
w("  results/statistical_outputs/data_quality_summary.txt")
w("")

cat("Inspection complete.\n")
cat("Project root: ", getwd(), "\n", sep = "")
cat("Wrote: results/tables/data_structure_summary.csv\n")
cat("Wrote: results/tables/missing_values_summary.csv\n")
cat("Wrote: results/statistical_outputs/data_quality_summary.txt\n")
