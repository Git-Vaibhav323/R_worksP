# Statistical Analysis of Rainfall Variability and Regional Monsoon Patterns in India

Undergraduate Probability and Statistics project in R. Phase 1 prepares a reproducible data-quality audit and cleaned analysis copies. Original government CSV files are never modified.

## Datasets (already in `data/`)

Both files are India Meteorological Department (IMD) rainfall series in **millimetres**. Zero means no rain that month; `NA` means the value was not reported.

| File | Spatial support | Time span in file | Identifying columns | Rainfall columns |
| --- | --- | --- | --- | --- |
| `data/rainfall_area-wt_India_1901-2015.csv` | One all-India **area-weighted** series (`REGION = INDIA`) | 1901–2015 | `REGION`, `YEAR` | Months `JAN`–`DEC`; seasons `Jan-Feb`, `Mar-May`, `Jun-Sep`, `Oct-Dec`; `ANNUAL` |
| `data/Sub_Division_IMD_2017.csv` | IMD **meteorological subdivisions** | 1901–2017 (unbalanced: not every subdivision has every year) | `SUBDIVISION`, `YEAR` | Months `JAN`–`DEC`; seasons `JF`, `MAM`, `JJAS`, `OND`; `ANNUAL` |

Season groups follow the standard IMD calendar:

- Winter: January–February
- Pre-monsoon: March–May
- Southwest monsoon: June–September
- Post-monsoon: October–December

The two files do **not** share the same spatial support. Subdivision rows must not be averaged with equal weights to reconstruct the national series.

Do not replace these files with Kaggle copies or other downloads.

## Project structure

```
data/                              original CSVs (read-only)
R_code/                            R scripts
  01_data_inspection.R
  02_data_cleaning.R
graphs/                            reserved for later figures
results/tables/                    CSV summaries and the cleaning log
results/statistical_outputs/       text audit
results/cleaned_data/              analysis copies only
report/                            reserved for the academic write-up
```

## How to run (Windows, RStudio)

1. Install R (this project was executed with R 4.6.1). RStudio is recommended.
2. Open the project folder as the working directory:
   - **Session → Set Working Directory → Choose Directory…** and select the folder that contains `data/` and `R_code/`, **or**
   - **File → Open Project** if you later add an `.Rproj` file.
3. Scripts use **relative paths** and also reset the working directory to the project root when run with `Rscript` from `R_code/`.
4. In the Console:

```r
source("R_code/01_data_inspection.R")
source("R_code/02_data_cleaning.R")
source("R_code/03_descriptive_analysis.R")
source("R_code/04_probability_distributions.R")
source("R_code/05_sampling_and_intervals.R")
source("R_code/06_hypothesis_testing.R")
source("R_code/07_correlation_regression.R")
source("R_code/08_regional_analysis.R")
source("R_code/09_time_series_analysis.R")
```

Or run Phase 2 in one call: `source("R_code/10_run_phase2.R")`.

From PowerShell in the project root (adjust the R path if needed):

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" R_code/01_data_inspection.R
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" R_code/02_data_cleaning.R
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" R_code/10_run_phase2.R
```

Base R only: no extra packages are required.

## Phase 1 outputs

- `results/statistical_outputs/data_quality_summary.txt`
- `results/tables/data_structure_summary.csv`
- `results/tables/missing_values_summary.csv`
- `results/tables/cleaning_log.csv`
- `results/tables/seasonal_annual_consistency.csv`
- `results/tables/subdivision_year_coverage.csv`
- `results/tables/subdivision_year_completeness.csv`
- `results/cleaned_data/india_national_rainfall_cleaned.csv`
- `results/cleaned_data/imd_subdivision_rainfall_cleaned.csv`

## Cleaning policy (Phase 1)

- Original files in `data/` are left unchanged.
- Incomplete years are **kept**; missing months stay `NA`. No listwise deletion of incomplete rows.
- Zeros are **not** treated as missing.
- No arbitrary imputation of unobserved monthly rainfall.
- Where a source seasonal or annual total is missing **and** all constituent months are present, the `*_analysis` columns contain the month sum. That reconstructs a derived total; it does not invent rain.
- Source seasonal/annual values are otherwise retained. Differences greater than 1 mm versus the month sum are **flagged**.
- Outliers are flagged with a conservative 3×IQR rule (by month nationally; by subdivision and month for the panel). Extreme years are **not** deleted.
- Subdivision names are trimmed; historical IMD labels such as Orissa, Matathwada and Rayalseema are kept.

## Known limitations (do not ignore in later phases)

- The subdivision table is an **unbalanced panel** (late starts, internal missing years such as Andaman & Nicobar Islands 1909).
- Many subdivision years have one or more monthly `NA`s; annual and monsoon analyses should use complete-month rows or the documented `*_analysis` fields.
- National series ends in 2015; subdivision series ends in 2017.
- Area-weighted India is not comparable to a simple mean of subdivisions.

Phase 2 scripts (`R_code/03`–`09`) produce descriptive statistics, probability summaries, distribution fits, sampling demonstrations, hypothesis tests, correlation, regression, regional ANOVA, and time-series trend diagnostics. See `results/tables/concept_coverage_matrix.csv`.

The academic report is `report/project_report.Rmd`. Knit it after the analysis outputs exist:

```r
source("R_code/11_render_report.R")
```

Deliverables: `report/Rainfall_Statistics_India_Project_Report.docx` and `.pdf`. See `report/README.md`.
