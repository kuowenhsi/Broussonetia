local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("vcfR", "ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "purrr", "forcats")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

library(vcfR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(scales)
library(purrr)
library(forcats)

vcf_file <- "VCF/filtered/Broussonetia_RADseq.09_missing50_mac1_site_missing20.vcf.gz"
metadata_file <- "reports/missing50_mac1/PRJNA437223_Broussonetia_meta.missing50_mac1.csv"
out_dir <- "reports/missing20_100kb_windows_ld"
window_size <- 100000L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading VCF: ", vcf_file)
vcf <- read.vcfR(vcf_file, verbose = FALSE)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
gt <- gsub("^([0-9.])\\|([0-9.])", "\\1/\\2", gt)
is_missing <- is.na(gt) | gt %in% c("./.", ".|.", ".")

index_stats <- system2("bcftools", c("index", "-s", vcf_file), stdout = TRUE) |>
  strsplit("\t") |>
  do.call(what = rbind) |>
  as_tibble(.name_repair = "minimal") |>
  setNames(c("chrom", "length", "sites")) |>
  mutate(
    length = as.integer(length),
    sites = as.integer(sites)
  )

chrom_info <- tibble(
  chrom = vcf@fix[, "CHROM"],
  pos = as.integer(vcf@fix[, "POS"])
) |>
  distinct(chrom, pos) |>
  group_by(chrom) |>
  summarise(max_pos = max(pos), .groups = "drop") |>
  mutate(chrom_num = as.integer(str_remove(chrom, "^scaffold_"))) |>
  filter(!is.na(chrom_num), chrom_num <= 18) |>
  left_join(index_stats, by = "chrom") |>
  arrange(chrom_num) |>
  transmute(
    chrom,
    chrom_num,
    chrom_label = paste0("Chr", chrom_num),
    length
  ) |>
  mutate(
    offset = lag(cumsum(length), default = 0),
    chrom_mid = offset + length / 2
  )

variants <- tibble(
  row_id = seq_len(nrow(vcf@fix)),
  chrom = vcf@fix[, "CHROM"],
  pos = as.integer(vcf@fix[, "POS"])
) |>
  inner_join(chrom_info, by = "chrom") |>
  mutate(
    window_start = ((pos - 1L) %/% window_size) * window_size + 1L,
    window_end = pmin(window_start + window_size - 1L, length),
    window_mid = (window_start + window_end) / 2,
    genome_window_mid = offset + window_mid
  )

chrom_windows <- chrom_info |>
  transmute(
    chrom,
    chrom_num,
    chrom_label,
    length,
    offset,
    starts = map(length, ~ seq(1L, .x, by = window_size))
  ) |>
  unnest(starts) |>
  rename(window_start = starts) |>
  mutate(
    window_end = pmin(window_start + window_size - 1L, length),
    window_mid = (window_start + window_end) / 2,
    genome_window_mid = offset + window_mid,
    window_width = window_end - window_start + 1L
  )

site_counts <- variants |>
  count(chrom, window_start, name = "site_count") |>
  right_join(chrom_windows, by = c("chrom", "window_start")) |>
  mutate(
    site_count = replace_na(site_count, 0L),
    sites_per_100kb = site_count / (window_width / window_size)
  ) |>
  arrange(chrom_num, window_start)

write_csv(site_counts, file.path(out_dir, "site_counts_100kb_continuous.csv"))

metadata <- read_csv(metadata_file, show_col_types = FALSE) |>
  filter(Run %in% colnames(gt)) |>
  mutate(
    geo_loc_name_country = replace_na(geo_loc_name_country, "Unknown"),
    geo_loc_name_country = if_else(geo_loc_name_country == "", "Unknown", geo_loc_name_country)
  )

metadata <- metadata |>
  arrange(match(Run, colnames(gt)))

stopifnot(identical(metadata$Run, colnames(gt)))

country_levels <- metadata |>
  count(geo_loc_name_country, sort = TRUE) |>
  pull(geo_loc_name_country)

country_missing <- map_dfr(country_levels, function(country) {
  sample_cols <- metadata$Run[metadata$geo_loc_name_country == country]
  variant_missing_rate <- rowMeans(is_missing[, sample_cols, drop = FALSE])
  variants |>
    mutate(
      geo_loc_name_country = country,
      variant_missing_rate = variant_missing_rate
    ) |>
    group_by(geo_loc_name_country, chrom, window_start) |>
    summarise(
      sites = n(),
      missing_rate = mean(variant_missing_rate),
      .groups = "drop"
    )
})

country_windows <- crossing(
  geo_loc_name_country = country_levels,
  chrom_windows
) |>
  left_join(country_missing, by = c("geo_loc_name_country", "chrom", "window_start")) |>
  mutate(
    sites = replace_na(sites, 0L),
    geo_loc_name_country = factor(geo_loc_name_country, levels = rev(country_levels))
  ) |>
  arrange(geo_loc_name_country, chrom_num, window_start)

write_csv(country_windows, file.path(out_dir, "site_missing_rate_by_country_100kb_continuous.csv"))

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

chrom_axis <- scale_x_continuous(
  breaks = chrom_info$chrom_mid,
  labels = chrom_info$chrom_label,
  expand = expansion(mult = c(0.002, 0.002))
)

p_counts <- ggplot(site_counts, aes(genome_window_mid, site_count)) +
  geom_col(width = window_size * 0.9, fill = "#4269A8") +
  geom_vline(xintercept = chrom_info$offset[-1], color = "grey70", linewidth = 0.25) +
  chrom_axis +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Concatenated chromosome position",
    y = "Variant sites per 100 kb",
    title = "Site counts across 18 Broussonetia chromosomes",
    subtitle = "VCF: missing50_mac1_site_missing20"
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "site_counts_100kb_continuous.png"),
  p_counts,
  width = 14,
  height = 5.5,
  dpi = 300
)

p_missing <- ggplot(country_windows, aes(genome_window_mid, missing_rate)) +
  geom_col(aes(fill = missing_rate), width = window_size * 0.9, na.rm = TRUE) +
  geom_vline(xintercept = chrom_info$offset[-1], color = "grey78", linewidth = 0.2) +
  facet_grid(rows = vars(geo_loc_name_country), switch = "y") +
  chrom_axis +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_viridis_c(option = "magma", labels = percent, limits = c(0, 1), na.value = "grey95") +
  labs(
    x = "Concatenated chromosome position",
    y = NULL,
    fill = "Missing",
    title = "Windowed site missingness by country",
    subtitle = "Mean per-site missing rate among samples from each country; 100 kb windows"
  ) +
  plot_theme +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8),
    panel.spacing.y = unit(0.08, "lines"),
    legend.position = "right"
  )

ggsave(
  file.path(out_dir, "site_missing_rate_by_country_100kb_continuous.png"),
  p_missing,
  width = 14,
  height = 12,
  dpi = 300
)

write_csv(
  chrom_info,
  file.path(out_dir, "chromosome_offsets_for_continuous_axis.csv")
)

message("Wrote window summaries and plots to: ", normalizePath(out_dir))
