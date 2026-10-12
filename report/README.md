# How to regenerate the academic report

The report is `report/project_report.Rmd`. Numbers are read from `results/tables/` at knit time, and figures are taken from `graphs/`. Do not type results by hand.

## Prerequisites

- R 4.x and RStudio
- R package `rmarkdown`
- Pandoc (bundled with RStudio). The render script looks for  
  `C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc.exe`
- Phase 1 and Phase 2 outputs already in `results/` and `graphs/`

PDF: `pdflatex` is not required. If LaTeX is missing, the render script asks Microsoft Word to save the DOCX as PDF (verified on this machine). You can also open the DOCX and use **File → Save As → PDF**.

## From the project root in RStudio

```r
# Optional: rebuild tables and graphs
source("R_code/01_data_inspection.R")
source("R_code/02_data_cleaning.R")
source("R_code/10_run_phase2.R")

# Knit DOCX and PDF
source("R_code/11_render_report.R")
```

From PowerShell:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" R_code/11_render_report.R
```

## Outputs

- `report/Rainfall_Statistics_India_Project_Report.docx`
- `report/Rainfall_Statistics_India_Project_Report.pdf` (Word Save-As PDF if `pdflatex` is not installed)
- `report/render_log.txt`

The Rmd lives in `report/` and loads figures from `../graphs/`. If you Knit from the RStudio button, keep the knit directory as the `report/` folder (the default for this file). To rebuild from the project root, use `source("R_code/11_render_report.R")`.

## Student details

Put name, roll number and institution on the title page before submission. Do not change plotted data or invented statistics.
