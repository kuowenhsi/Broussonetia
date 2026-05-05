local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("vcfR", "ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "purrr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

library(vcfR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(scales)

vcf_file <- "VCF/Broussonetia_RADseq.07_tags_refilled.vcf.gz"
fai_file <- "genome_AS/B.papyrifera_yahs.out_scaffolds.fa.fai"
out_dir <- "reports/variant_distribution_scaffolds"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

window_sizes <- c(100000, 250000, 500000, 1000000)
recommended_window <- 500000

chrom_info <- read_tsv(
  fai_file,
  col_names = c("chrom", "length", "offset", "line_bases", "line_width"),
  show_col_types = FALSE
) |>
  slice_head(n = 18) |>
  mutate(
    chrom_num = as.integer(str_remove(chrom, "^scaffold_")),
    chrom = factor(chrom, levels = paste0("scaffold_", 1:18)),
    chrom_label = paste0("Chr", chrom_num),
    chrom_label = factor(chrom_label, levels = paste0("Chr", 1:18)),
    length_mb = length / 1e6
  )

message("Reading VCF with vcfR: ", vcf_file)
vcf <- read.vcfR(vcf_file, verbose = FALSE)

variants <- tibble(
  chrom = vcf@fix[, "CHROM"],
  pos = as.integer(vcf@fix[, "POS"])
) |>
  filter(chrom %in% as.character(chrom_info$chrom)) |>
  mutate(
    chrom = factor(chrom, levels = levels(chrom_info$chrom)),
    chrom_label = factor(
      paste0("Chr", str_remove(as.character(chrom), "^scaffold_")),
      levels = levels(chrom_info$chrom_label)
    )
  )

variant_counts <- chrom_info |>
  count(chrom, chrom_num, chrom_label, length, length_mb, name = "dummy") |>
  select(-dummy) |>
  left_join(
    variants |> count(chrom, name = "variant_count"),
    by = "chrom"
  ) |>
  mutate(
    variant_count = replace_na(variant_count, 0L),
    variants_per_mb = variant_count / length_mb
  )

write_csv(variant_counts, file.path(out_dir, "variant_counts_by_chromosome.csv"))

make_window_counts <- function(window_size) {
  chrom_windows <- chrom_info |>
    transmute(
      chrom,
      chrom_label,
      chrom_length = length,
      window_start = purrr::map(length, ~ seq(1L, .x, by = window_size))
    ) |>
    unnest(window_start) |>
    mutate(
      window_end = pmin(window_start + window_size - 1L, chrom_length),
      window_mid = (window_start + window_end) / 2,
      window_width = window_end - window_start + 1L,
      window_size = window_size
    )

  variants |>
    mutate(window_start = ((pos - 1L) %/% window_size) * window_size + 1L) |>
    count(chrom, window_start, name = "variant_count") |>
    right_join(chrom_windows, by = c("chrom", "window_start")) |>
    mutate(
      variant_count = replace_na(variant_count, 0L),
      variants_per_mb = variant_count / (window_width / 1e6),
      window_size_label = case_when(
        window_size >= 1e6 ~ paste0(window_size / 1e6, " Mb"),
        TRUE ~ paste0(window_size / 1000, " kb")
      ),
      window_size_label = factor(
        window_size_label,
        levels = c("100 kb", "250 kb", "500 kb", "1 Mb")
      )
    )
}

window_counts <- bind_rows(lapply(window_sizes, make_window_counts))
write_csv(window_counts, file.path(out_dir, "variant_counts_by_window_all_sizes.csv"))

best_counts <- window_counts |>
  filter(window_size == recommended_window)
write_csv(
  best_counts,
  file.path(out_dir, paste0("variant_counts_by_window_", recommended_window / 1000, "kb.csv"))
)

plot_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    axis.text.y = element_text(size = 9)
  )

p_counts <- variant_counts |>
  ggplot(aes(x = reorder(chrom_label, chrom_num), y = variant_count)) +
  geom_col(fill = "#4C78A8", width = 0.75) +
  geom_text(aes(label = comma(variant_count)), hjust = -0.08, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL,
    y = "Variant sites",
    title = "Variant sites per chromosome-scale scaffold",
    subtitle = "First 18 scaffolds from the FASTA index are treated as chromosomes"
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "variant_counts_by_chromosome.png"),
  p_counts,
  width = 8,
  height = 6,
  dpi = 300
)

p_window_compare <- window_counts |>
  ggplot(aes(x = window_mid / 1e6, y = chrom_label, fill = variants_per_mb)) +
  geom_tile(aes(width = window_width / 1e6), height = 0.82) +
  facet_wrap(~ window_size_label, ncol = 1) +
  scale_fill_viridis_c(option = "magma", labels = comma, trans = "sqrt") +
  labs(
    x = "Position (Mb)",
    y = NULL,
    fill = "Variants / Mb",
    title = "Variant-site density across chromosomes at four window sizes",
    subtitle = "Smaller windows show local RAD tag clustering; larger windows smooth the sparse chromosomes"
  ) +
  plot_theme +
  theme(legend.position = "right")

ggsave(
  file.path(out_dir, "variant_density_window_size_comparison.png"),
  p_window_compare,
  width = 11,
  height = 11,
  dpi = 300
)

p_recommended <- best_counts |>
  ggplot(aes(x = window_mid / 1e6, y = chrom_label, fill = variants_per_mb)) +
  geom_tile(aes(width = window_width / 1e6), height = 0.82) +
  scale_fill_viridis_c(option = "magma", labels = comma, trans = "sqrt") +
  labs(
    x = "Position (Mb)",
    y = NULL,
    fill = "Variants / Mb",
    title = "Variant-site density across 18 Broussonetia chromosomes",
    subtitle = paste0("Recommended exploratory window: ", recommended_window / 1000, " kb")
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "variant_density_recommended_500kb.png"),
  p_recommended,
  width = 11,
  height = 6.5,
  dpi = 300
)

p_lines <- best_counts |>
  ggplot(aes(x = window_mid / 1e6, y = variant_count)) +
  geom_col(fill = "#4C78A8", width = 0.45) +
  facet_wrap(~ chrom_label, scales = "free_x", ncol = 3) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Position (Mb)",
    y = "Variant sites per 500 kb",
    title = "Per-chromosome variant-site distribution",
    subtitle = "Faceted view using the recommended 500 kb window"
  ) +
  plot_theme

ggsave(
  file.path(out_dir, "variant_distribution_faceted_500kb.png"),
  p_lines,
  width = 12,
  height = 10,
  dpi = 300
)

window_summary <- window_counts |>
  group_by(window_size, window_size_label) |>
  summarise(
    windows = n(),
    median_variants = median(variant_count),
    mean_variants = mean(variant_count),
    max_variants = max(variant_count),
    empty_windows = sum(variant_count == 0),
    empty_window_fraction = empty_windows / windows,
    .groups = "drop"
  )

write_csv(window_summary, file.path(out_dir, "window_size_summary.csv"))

message("Wrote outputs to: ", normalizePath(out_dir))
message("Recommended first visual window: ", recommended_window / 1000, " kb")
