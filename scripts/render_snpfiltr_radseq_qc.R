local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c(
  "SNPfiltR", "vcfR", "tidyverse", "rmarkdown", "knitr", "DT",
  "patchwork", "scales"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing, collapse = ", "),
    "\nInstall missing packages, then rerun this script.",
    call. = FALSE
  )
}

dir.create("reports", showWarnings = FALSE, recursive = TRUE)

rmarkdown::render(
  input = "reports/snpfiltr_radseq_qc.Rmd",
  output_file = "snpfiltr_radseq_qc.html",
  output_dir = "reports",
  clean = TRUE,
  envir = new.env(parent = globalenv())
)
