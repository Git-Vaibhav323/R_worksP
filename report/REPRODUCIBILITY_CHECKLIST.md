# Reproducibility checklist

Use this list before a viva or submission.

## Environment

- [ ] R version recorded (project executed with 4.6.1)
- [ ] Working directory is the project root (contains `data/` and `R_code/`)
- [ ] No extra rainfall files were added to `data/`
- [ ] Original CSVs in `data/` were not edited

## Phase 1

- [ ] `source("R_code/01_data_inspection.R")` completes
- [ ] `source("R_code/02_data_cleaning.R")` completes
- [ ] `results/cleaned_data/india_national_rainfall_cleaned.csv` exists (115 rows)
- [ ] `results/cleaned_data/imd_subdivision_rainfall_cleaned.csv` exists (4188 rows)
- [ ] `results/tables/cleaning_log.csv` exists

## Phase 2

- [ ] `source("R_code/10_run_phase2.R")` completes (or `03`–`09` in order)
- [ ] Random seed 19012015 used for CLT and bootstrap KS
- [ ] Twenty PNG files exist in `graphs/`
- [ ] Concept matrix exists: `results/tables/concept_coverage_matrix.csv`

## Report

- [ ] `source("R_code/11_render_report.R")` writes the DOCX
- [ ] Figures in the Word file match `graphs/` (not screenshots from the web)
- [ ] Inline numbers match `results/tables/` (the Rmd reads those CSVs)
- [ ] Student name and institution filled on the title page
- [ ] PDF produced, **or** Word used to export PDF, with the limitation stated

## Scientific honesty

- [ ] Zeros are not treated as missing
- [ ] Incomplete subdivision years are not secretly dropped from every analysis (annual regional work uses complete months)
- [ ] National area-weighted series is not replaced by a mean of subdivisions
- [ ] Correlation is not described as causation
- [ ] Non-rejection of H0 is not described as proof of H0
- [ ] Trend tests are not described as proof of climate change
- [ ] No Kaggle or invented observations
