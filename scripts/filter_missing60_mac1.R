local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("vcfR", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

library(vcfR)
library(data.table)

original_vcf <- "VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz"
metadata_file <- "PRJNA437223_Broussonetia_meta.csv"
out_vcf <- "VCF/filtered/Broussonetia_RADseq.08_missing60_mac1.vcf.gz"
report_dir <- "reports/missing60_mac1"
sample_qc_file <- file.path(report_dir, "sample_missingness_before_missing60_filter.csv")
retained_samples_file <- file.path(report_dir, "samples_retained_missing60.txt")
removed_samples_file <- file.path(report_dir, "samples_removed_missing60.txt")
filtered_metadata_file <- file.path(report_dir, "PRJNA437223_Broussonetia_meta.missing60_mac1.csv")
summary_file <- file.path(report_dir, "missing60_mac1_filter_summary.csv")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_vcf), recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
}

stop_if_missing(original_vcf)
stop_if_missing(metadata_file)

needs_refresh <- function(target, sources) {
  if (!file.exists(target)) {
    return(TRUE)
  }
  target_mtime <- file.info(target)$mtime
  any(vapply(sources, function(source) {
    file.exists(source) && file.info(source)$mtime > target_mtime
  }, logical(1)))
}

if (needs_refresh(sample_qc_file, original_vcf) ||
    needs_refresh(retained_samples_file, c(original_vcf, sample_qc_file)) ||
    needs_refresh(removed_samples_file, c(original_vcf, sample_qc_file))) {
  message("Reading source VCF to calculate sample missingness: ", original_vcf)
  v <- read.vcfR(original_vcf, verbose = FALSE)
  gt <- extract.gt(v, element = "GT", as.numeric = FALSE)
  called <- !is.na(gt) & !(gt %in% c("./.", ".|.", "."))
  missing_rate <- colMeans(!called)
  sample_qc <- data.table(
    Run = names(missing_rate),
    missing_rate = as.numeric(missing_rate),
    retain_missing60 = missing_rate <= 0.60
  )
  fwrite(sample_qc, sample_qc_file)
  writeLines(sample_qc[retain_missing60 == TRUE, Run], retained_samples_file)
  writeLines(sample_qc[retain_missing60 == FALSE, Run], removed_samples_file)
}

if (needs_refresh(out_vcf, c(original_vcf, retained_samples_file))) {
  message("Filtering retained samples and MAC >= 1 sites")
  cmd <- paste(
    "bcftools view -S", shQuote(retained_samples_file), "-Ou", shQuote(original_vcf),
    "| bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,NS,F_MISSING,HWE,ExcHet",
    "| bcftools view -i 'MAF>0' -Oz -o", shQuote(out_vcf)
  )
  if (system(cmd) != 0) {
    stop("Could not generate filtered VCF: ", out_vcf, call. = FALSE)
  }
  if (system2("bcftools", c("index", "-f", out_vcf)) != 0) {
    stop("Could not index filtered VCF: ", out_vcf, call. = FALSE)
  }
}

retained_samples <- fread(retained_samples_file, header = FALSE)$V1
metadata <- fread(metadata_file)
fwrite(metadata[Run %in% retained_samples], filtered_metadata_file)

source_sites <- as.integer(system2("bcftools", c("index", "-n", original_vcf), stdout = TRUE))
filtered_sites <- as.integer(system2("bcftools", c("index", "-n", out_vcf), stdout = TRUE))
summary_dt <- data.table(
  source_vcf = original_vcf,
  filtered_vcf = out_vcf,
  source_sites = source_sites,
  filtered_sites_mac1 = filtered_sites,
  source_samples = length(system2("bcftools", c("query", "-l", original_vcf), stdout = TRUE)),
  retained_samples = length(retained_samples),
  removed_samples = length(fread(removed_samples_file, header = FALSE)$V1),
  sample_missing_filter = "remove samples with missing_rate > 0.60",
  site_filter = "MAF > 0 after retaining samples; equivalent to MAC >= 1 for polymorphic sites"
)
fwrite(summary_dt, summary_file)

message("Done. Output VCF: ", out_vcf)
