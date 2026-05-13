local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("vcfR", "tidyverse", "rmarkdown", "knitr", "DT", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

dir.create("reports/missing50_mac1", showWarnings = FALSE, recursive = TRUE)

project_dir <- getwd()
original_vcf <- file.path(project_dir, "VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz")
filtered_vcf_project <- file.path(project_dir, "VCF/filtered/Broussonetia_RADseq.08_missing50_mac1.vcf.gz")
site_missing25_vcf <- file.path(project_dir, "VCF/filtered/Broussonetia_RADseq.09_missing50_mac1_site_missing25.vcf.gz")
site_missing20_vcf <- file.path(project_dir, "VCF/filtered/Broussonetia_RADseq.09_missing50_mac1_site_missing20.vcf.gz")
filtered_vcf_temp <- "/private/tmp/Broussonetia_RADseq.08_missing50_mac1.vcf.gz"
sample_qc_file <- file.path(project_dir, "reports/missing50_mac1/sample_missingness_before_missing50_filter.csv")
retained_samples_file <- file.path(project_dir, "reports/missing50_mac1/samples_retained_missing50.txt")
removed_samples_file <- file.path(project_dir, "reports/missing50_mac1/samples_removed_missing50.txt")

if (!file.exists(original_vcf)) {
  stop("Original VCF does not exist: ", original_vcf, call. = FALSE)
}

needs_refresh <- function(target, sources) {
  if (!file.exists(target)) {
    return(TRUE)
  }
  target_mtime <- file.info(target)$mtime
  any(vapply(sources, function(source) {
    file.exists(source) && file.info(source)$mtime > target_mtime
  }, logical(1)))
}

if (needs_refresh(retained_samples_file, c(original_vcf, sample_qc_file)) ||
    needs_refresh(sample_qc_file, original_vcf) ||
    needs_refresh(removed_samples_file, c(original_vcf, sample_qc_file))) {
  v <- vcfR::read.vcfR(original_vcf, verbose = FALSE)
  gt <- vcfR::extract.gt(v, element = "GT", as.numeric = FALSE)
  called <- !is.na(gt) & gt != "./." & gt != ".|." & gt != "."
  missing_rate <- colMeans(!called)
  sample_qc <- data.frame(
    Run = names(missing_rate),
    missing_rate = as.numeric(missing_rate),
    retain_missing50 = missing_rate <= 0.5
  )
  write.csv(sample_qc, sample_qc_file, row.names = FALSE)
  writeLines(sample_qc$Run[sample_qc$retain_missing50], retained_samples_file)
  writeLines(sample_qc$Run[!sample_qc$retain_missing50], removed_samples_file)
}

if (needs_refresh(filtered_vcf_project, c(original_vcf, retained_samples_file))) {
  dir.create(dirname(filtered_vcf_project), showWarnings = FALSE, recursive = TRUE)
  cmd <- paste(
    "bcftools view -S", shQuote(retained_samples_file), "-Ou", shQuote(original_vcf),
    "| bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,NS,HWE,ExcHet",
    "| bcftools view -i 'MAF>0' -Oz -o", shQuote(filtered_vcf_project)
  )
  if (system(cmd) != 0) {
    stop("Could not regenerate filtered VCF with sample missingness filter.", call. = FALSE)
  }
  if (system2("bcftools", c("index", "-f", filtered_vcf_project)) != 0) {
    stop("Could not index filtered VCF: ", filtered_vcf_project, call. = FALSE)
  }
}

filter_site_missingness <- function(source_vcf, target_vcf, max_missing) {
  if (!needs_refresh(target_vcf, c(source_vcf, paste0(source_vcf, ".csi")))) {
    return(invisible(FALSE))
  }
  dir.create(dirname(target_vcf), showWarnings = FALSE, recursive = TRUE)
  expr <- sprintf("F_MISSING<%s && MAF>0", max_missing)
  if (system2("bcftools", c("view", "-i", shQuote(expr), "-Oz", "-o", target_vcf, source_vcf)) != 0) {
    stop("Could not filter sites with expression: ", expr, call. = FALSE)
  }
  if (system2("bcftools", c("index", "-f", target_vcf)) != 0) {
    stop("Could not index site-filtered VCF: ", target_vcf, call. = FALSE)
  }
  invisible(TRUE)
}

filter_site_missingness(filtered_vcf_project, site_missing25_vcf, 0.25)
filter_site_missingness(filtered_vcf_project, site_missing20_vcf, 0.20)

copy_if_needed <- function(from, to) {
  if (!file.exists(from)) {
    return(invisible(FALSE))
  }
  if (!file.exists(to) || file.info(from)$mtime > file.info(to)$mtime) {
    ok <- file.copy(from, to, overwrite = TRUE)
    if (!ok) {
      stop("Could not copy ", from, " to ", to, call. = FALSE)
    }
  }
  invisible(TRUE)
}

copy_if_needed(filtered_vcf_project, filtered_vcf_temp)
copy_if_needed(paste0(filtered_vcf_project, ".csi"), paste0(filtered_vcf_temp, ".csi"))

rmarkdown::render(
  input = "reports/missing50_mac1/missing50_mac1_hwe_report.Rmd",
  output_file = "missing50_mac1_hwe_report.html",
  output_dir = "reports/missing50_mac1",
  params = list(
    original_vcf = file.path(project_dir, "VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz"),
    filtered_vcf = filtered_vcf_temp,
    filtered_vcf_project = filtered_vcf_project,
    site_missing25_vcf = site_missing25_vcf,
    site_missing20_vcf = site_missing20_vcf,
    metadata_file = file.path(project_dir, "PRJNA437223_Broussonetia_meta.csv"),
    filtered_metadata_out = file.path(project_dir, "reports/missing50_mac1/PRJNA437223_Broussonetia_meta.missing50_mac1.csv")
  ),
  clean = TRUE,
  envir = new.env(parent = globalenv())
)
