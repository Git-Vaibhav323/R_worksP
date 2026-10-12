# Render the academic report to DOCX (and PDF if a LaTeX engine is available).
# Run from the project root after Phase 1–2 scripts have produced tables and graphs.

ensure <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
    d <- dirname(script_path)
    if (basename(d) == "R_code") setwd(dirname(d))
  } else if (basename(getwd()) == "R_code") setwd("..")
}
ensure()

pandoc_dir <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
if (file.exists(file.path(pandoc_dir, "pandoc.exe"))) {
  Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required. Install it in RStudio with install.packages('rmarkdown').")
}
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc was not found. Set RSTUDIO_PANDOC to the folder that contains pandoc.exe.")
}

rmd <- file.path("report", "project_report.Rmd")
if (!file.exists(rmd)) stop("Missing ", rmd)

message("Pandoc: ", rmarkdown::pandoc_version())
message("Working directory: ", getwd())

root <- normalizePath(getwd(), winslash = "/")
report_dir <- file.path(root, "report")
resource_path <- paste(c(report_dir, root, file.path(root, "graphs")), collapse = ";")
word_fmt <- rmarkdown::word_document(
  toc = TRUE,
  toc_depth = 3,
  fig_caption = TRUE,
  pandoc_args = c(paste0("--resource-path=", resource_path))
)

# Knit from report/ so ../graphs and ../results resolve for both knitr and pandoc.
docx <- rmarkdown::render(
  input = rmd,
  output_format = word_fmt,
  output_file = "Rainfall_Statistics_India_Project_Report.docx",
  output_dir = report_dir,
  knit_root_dir = report_dir,
  quiet = FALSE,
  envir = new.env()
)
message("Wrote: ", docx)
if (file.exists(docx) && file.info(docx)$size < 200000) {
  warning("DOCX is unusually small; figures may be missing.")
}

pdf_ok <- FALSE
pdf_error <- NULL
pdf_path <- file.path("report", "Rainfall_Statistics_India_Project_Report.pdf")

tryCatch({
  pdf <- rmarkdown::render(
    input = rmd,
    output_format = rmarkdown::pdf_document(
      toc = TRUE,
      toc_depth = 3,
      number_sections = TRUE,
      fig_caption = TRUE,
      pandoc_args = c(paste0("--resource-path=", resource_path))
    ),
    output_file = "Rainfall_Statistics_India_Project_Report.pdf",
    output_dir = report_dir,
    knit_root_dir = report_dir,
    quiet = FALSE,
    envir = new.env()
  )
  pdf_ok <- file.exists(pdf)
  if (pdf_ok) message("Wrote: ", pdf)
}, error = function(e) {
  pdf_error <<- conditionMessage(e)
  message("LaTeX PDF skipped: ", pdf_error)
})

# Fallback: Microsoft Word Save-As PDF if pdflatex is absent.
if (!pdf_ok && file.exists(docx)) {
  word_exe <- "C:/Program Files/Microsoft Office/root/Office16/WINWORD.EXE"
  ps <- paste0(
    "$docx = (Resolve-Path '", normalizePath(docx, winslash = "/"), "').Path; ",
    "$pdf = [System.IO.Path]::ChangeExtension($docx, '.pdf'); ",
    "$w = $null; try { $w = New-Object -ComObject Word.Application; ",
    "$w.Visible = $false; $d = $w.Documents.Open($docx); ",
    "$d.SaveAs([ref]$pdf, [ref]17); $d.Close(); $w.Quit(); ",
    "Write-Output $pdf } catch { if ($w) { $w.Quit() }; throw }"
  )
  tryCatch({
    out <- system2("powershell", c("-NoProfile", "-Command", ps), stdout = TRUE, stderr = TRUE)
    if (file.exists(pdf_path)) {
      pdf_ok <- TRUE
      pdf_error <- paste("PDF produced via Word COM, not pdflatex.", paste(out, collapse = " "))
      message("Wrote PDF via Word: ", pdf_path)
    } else {
      pdf_error <- paste(pdf_error, "Word COM did not create PDF.", paste(out, collapse = " "))
    }
  }, error = function(e) {
    pdf_error <<- paste(pdf_error, "Word COM failed:", conditionMessage(e))
  })
}

writeLines(
  c(
    paste("Report render:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("DOCX:", docx, "exists=", file.exists(docx)),
    paste("PDF attempted:", TRUE),
    paste("PDF produced:", pdf_ok),
    paste("PDF error:", ifelse(is.null(pdf_error), "none", pdf_error))
  ),
  file.path("report", "render_log.txt")
)
