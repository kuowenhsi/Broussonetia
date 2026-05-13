local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("vcfR", "ape", "ggplot2", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

library(vcfR)
library(ape)
library(ggplot2)
library(data.table)

vcf_in <- "VCF/filtered/Broussonetia_RADseq.08_missing60_mac1.vcf.gz"
metadata_file <- "reports/missing60_mac1/PRJNA437223_Broussonetia_meta.missing60_mac1.csv"
out_dir <- "reports/missing60_mac1_site_missing_pca_nj"
vcf_dir <- file.path(out_dir, "vcf")
thresholds <- c(80, 60, 40, 20)
outgroup_taxa <- c("Broussonetia kaempferi", "Broussonetia monoica")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vcf_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
}

stop_if_missing(vcf_in)
stop_if_missing(metadata_file)

message("Reading metadata: ", metadata_file)
metadata <- fread(metadata_file)
needed_cols <- c("Run", "Organism", "geo_loc_name_country", "geo_loc_name", "haplotype")
if (!all(needed_cols %in% names(metadata))) {
  stop("Metadata must contain columns: ", paste(needed_cols, collapse = ", "), call. = FALSE)
}

clean_location <- function(x) {
  x <- gsub("\\\\,", ",", x)
  x <- gsub(":\\s*", ", ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

metadata[, country := fifelse(
  is.na(geo_loc_name_country) | geo_loc_name_country == "",
  "Unknown",
  geo_loc_name_country
)]
metadata[, locality := fifelse(
  is.na(geo_loc_name) | geo_loc_name == "",
  country,
  clean_location(geo_loc_name)
)]
metadata[, haplotype_clean := fifelse(
  is.na(haplotype) | haplotype == "" | haplotype == "missing",
  NA_character_,
  as.character(haplotype)
)]
metadata[, tree_label := fifelse(
  is.na(haplotype_clean),
  locality,
  paste0(locality, " [", haplotype_clean, "]")
)]

samples <- system2("bcftools", c("query", "-l", vcf_in), stdout = TRUE)
metadata <- metadata[Run %in% samples]
metadata <- metadata[match(samples, metadata$Run)]

if (!identical(metadata$Run, samples)) {
  missing_meta <- samples[!samples %in% metadata$Run]
  stop("Missing metadata for VCF sample(s): ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

outgroup_samples <- metadata[Organism %in% outgroup_taxa, Run]
if (length(outgroup_samples) != length(outgroup_taxa)) {
  found <- metadata[Organism %in% outgroup_taxa, .(Run, Organism)]
  stop(
    "Expected one retained sample for each outgroup taxon. Found: ",
    paste(found$Organism, found$Run, sep = "=", collapse = "; "),
    call. = FALSE
  )
}

country_levels <- metadata[, .N, by = country][order(-N, country)]$country
metadata[, country := factor(country, levels = country_levels)]
country_palette <- setNames(
  grDevices::hcl.colors(length(country_levels), palette = "Dark 3"),
  country_levels
)

gt_to_dosage <- function(gt) {
  gt <- gsub("\\|", "/", gt)
  dosage <- matrix(NA_real_, nrow = nrow(gt), ncol = ncol(gt), dimnames = dimnames(gt))
  dosage[gt == "0/0"] <- 0
  dosage[gt %in% c("0/1", "1/0")] <- 1
  dosage[gt == "1/1"] <- 2

  unresolved <- is.na(dosage) & !(is.na(gt) | gt %in% c("./.", ".", ".|."))
  if (any(unresolved)) {
    values <- gt[unresolved]
    generic <- vapply(strsplit(values, "/", fixed = TRUE), function(a) {
      if (length(a) == 0 || any(a == ".")) {
        return(NA_real_)
      }
      sum(suppressWarnings(as.integer(a)) > 0, na.rm = TRUE)
    }, numeric(1))
    dosage[unresolved] <- generic
  }

  storage.mode(dosage) <- "double"
  dosage
}

impute_site_means <- function(mat) {
  site_means <- rowMeans(mat, na.rm = TRUE)
  site_means[!is.finite(site_means)] <- 0
  miss <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(miss) > 0) {
    mat[miss] <- site_means[miss[, "row"]]
  }
  mat
}

pairwise_genetic_distance <- function(mat) {
  present <- !is.na(mat)
  x <- mat
  x[!present] <- 0
  m <- matrix(as.numeric(present), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

  shared_sites <- crossprod(m)
  cross <- crossprod(x)
  ss_left <- crossprod(x^2, m)
  ss_right <- t(ss_left)
  d2 <- (ss_left + ss_right - 2 * cross) / shared_sites
  d2[shared_sites == 0] <- NA_real_
  d2[d2 < 0 & d2 > -1e-10] <- 0
  d <- sqrt(d2)
  diag(d) <- 0
  as.dist(d)
}

root_on_outgroups <- function(tree) {
  root(
    tree,
    outgroup = outgroup_samples,
    resolve.root = TRUE
  )
}

plot_tree <- function(tree, tip_colors, threshold, path, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") {
    png(path, width = 6000, height = 7600, res = 300)
  } else {
    pdf(path, width = 20, height = 25)
  }
  on.exit(dev.off(), add = TRUE)
  par(mar = c(1, 1, 3, 1))
  plot(
    tree,
    type = "phylogram",
    direction = "rightwards",
    cex = 0.45,
    font = 1,
    tip.color = tip_colors,
    label.offset = 0.002,
    main = sprintf("Rooted neighbor-joining tree, site missing <= %d%%", threshold)
  )
  legend(
    "topleft",
    legend = names(country_palette),
    col = country_palette,
    pch = 19,
    cex = 0.65,
    bty = "n",
    title = "Country"
  )
}

run_threshold <- function(threshold) {
  tag <- paste0("site_missing", threshold)
  max_missing <- threshold / 100
  filtered_vcf <- file.path(vcf_dir, paste0("Broussonetia_RADseq.08_missing60_mac1_", tag, ".vcf.gz"))

  if (!file.exists(filtered_vcf)) {
    message("Filtering ", tag, " VCF")
    status <- system2(
      "bcftools",
      c(
        "view",
        "-i", shQuote(sprintf("F_MISSING<=%.2f", max_missing)),
        "-Oz",
        "-o", filtered_vcf,
        vcf_in
      )
    )
    if (!identical(status, 0L)) {
      stop("bcftools view failed for ", tag, call. = FALSE)
    }
  }

  index_file <- paste0(filtered_vcf, ".csi")
  if (!file.exists(index_file)) {
    status <- system2("bcftools", c("index", "-f", filtered_vcf))
    if (!identical(status, 0L)) {
      stop("bcftools index failed for ", tag, call. = FALSE)
    }
  }

  site_count <- as.integer(system2("bcftools", c("index", "-n", filtered_vcf), stdout = TRUE))
  message("Reading ", tag, ": ", site_count, " sites")
  vcf <- read.vcfR(filtered_vcf, verbose = FALSE)
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  gt <- gt[, samples, drop = FALSE]
  dosage <- gt_to_dosage(gt)

  missing_rate <- mean(is.na(dosage))
  pca_input <- impute_site_means(dosage)
  pca <- prcomp(t(pca_input), center = TRUE, scale. = FALSE)
  pve <- (pca$sdev^2) / sum(pca$sdev^2)

  pca_dt <- data.table(
    Run = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3]
  )
  pca_dt <- merge(
    pca_dt,
    metadata[, .(Run, Organism, country, locality, haplotype, tree_label)],
    by = "Run",
    sort = FALSE
  )
  fwrite(pca_dt, file.path(out_dir, paste0(tag, "_pca_coordinates.csv")))

  pca_plot <- ggplot(pca_dt, aes(PC1, PC2, color = country)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = country_palette, drop = FALSE) +
    labs(
      x = sprintf("PC1 (%.1f%%)", 100 * pve[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * pve[2]),
      color = "Country",
      title = sprintf("Broussonetia RADseq PCA, missing60/mac1, site missing <= %d%%", threshold),
      subtitle = sprintf("%s SNPs; %.2f%% genotype calls missing before mean imputation", format(site_count, big.mark = ","), 100 * missing_rate)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title.position = "plot"
    )

  ggsave(file.path(out_dir, paste0(tag, "_pca_country.png")), pca_plot, width = 8.5, height = 6.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(tag, "_pca_country.pdf")), pca_plot, width = 8.5, height = 6.5)

  message("Calculating rooted NJ tree for ", tag)
  dist_obj <- pairwise_genetic_distance(dosage)
  tree <- nj(dist_obj)
  tree <- root_on_outgroups(tree)
  tree_sample_labels <- tree$tip.label
  tree$tip.label <- metadata$tree_label[match(tree_sample_labels, metadata$Run)]
  write.tree(tree, file.path(out_dir, paste0(tag, "_rooted_nj_tree.newick")))

  tip_colors <- country_palette[as.character(metadata$country[match(tree_sample_labels, metadata$Run)])]
  tip_colors[is.na(tip_colors)] <- "grey30"
  names(tip_colors) <- tree$tip.label

  plot_tree(
    tree,
    tip_colors,
    threshold,
    file.path(out_dir, paste0(tag, "_rooted_nj_tree_location_labels.png")),
    "png"
  )
  plot_tree(
    tree,
    tip_colors,
    threshold,
    file.path(out_dir, paste0(tag, "_rooted_nj_tree_location_labels.pdf")),
    "pdf"
  )

  data.table(
    site_missing_filter_percent = threshold,
    sites_retained = site_count,
    samples = ncol(dosage),
    genotype_missing_rate = missing_rate,
    pc1_percent = 100 * pve[1],
    pc2_percent = 100 * pve[2],
    filtered_vcf = filtered_vcf,
    outgroup_samples = paste(outgroup_samples, collapse = ";")
  )
}

summary_dt <- rbindlist(lapply(thresholds, run_threshold))
fwrite(summary_dt, file.path(out_dir, "site_missing_filter_summary.csv"))

message("Done. Outputs written to ", out_dir)
