local_lib <- file.path(getwd(), "R_libs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c("ggplot2", "dplyr", "readr", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

library(ggplot2)
library(dplyr)
library(readr)
library(scales)

out_dir <- "reports/missing20_100kb_windows_ld/ld_decay_plink2"
ld_file <- file.path(out_dir, "missing20.ld_0_1000kb.vcor")
plink1_ld_file <- file.path(out_dir, "missing20.ld_0_1000kb_plink1.ld.gz")

if (file.exists(plink1_ld_file)) {
  ld_file <- plink1_ld_file
} else if (!file.exists(ld_file)) {
  candidates <- list.files(out_dir, pattern = "\\.vcor$", full.names = TRUE)
  if (length(candidates) == 1) {
    ld_file <- candidates
  } else {
    stop("Could not find PLINK2 LD output .vcor file in: ", out_dir, call. = FALSE)
  }
}

ld <- read_table(ld_file, show_col_types = FALSE)
names(ld) <- sub("^#", "", names(ld))

if (all(c("BP_A", "BP_B") %in% names(ld)) && !"DIST" %in% names(ld)) {
  ld <- ld |>
    mutate(DIST = abs(BP_B - BP_A))
}

dist_col <- intersect(c("DIST", "BP_DISTANCE", "DISTANCE"), names(ld))[1]
r2_col <- intersect(c("R2", "UNPHASED_R2"), names(ld))[1]

if (is.na(dist_col) || is.na(r2_col)) {
  stop(
    "Could not identify distance/R2 columns in ",
    ld_file,
    ". Columns are: ",
    paste(names(ld), collapse = ", "),
    call. = FALSE
  )
}

bin_size <- 1000L
max_dist <- 150000L

ld_decay <- ld |>
  transmute(
    distance_bp = as.numeric(.data[[dist_col]]),
    r2 = as.numeric(.data[[r2_col]])
  ) |>
  filter(!is.na(distance_bp), !is.na(r2), distance_bp >= 0, distance_bp <= max_dist) |>
  mutate(
    distance_bin_start = pmin(floor(distance_bp / bin_size) * bin_size, max_dist - bin_size),
    distance_bin_mid_kb = (distance_bin_start + bin_size / 2) / 1000
  ) |>
  group_by(distance_bin_start, distance_bin_mid_kb) |>
  summarise(
    pairs = n(),
    mean_r2 = mean(r2),
    median_r2 = median(r2),
    q25_r2 = quantile(r2, 0.25),
    q75_r2 = quantile(r2, 0.75),
    .groups = "drop"
  ) |>
  arrange(distance_bin_start)

write_csv(ld_decay, file.path(out_dir, "ld_decay_0_1000kb_10kb_bins.csv"))

p <- ggplot(ld_decay, aes(x = distance_bin_mid_kb, y = mean_r2)) +
  # geom_point(data = ld, aes(x = DIST/1000, y = R2), color = "gray", alpha = 0.75) +
  geom_line(color = "red", linewidth = 0.8) +
  scale_x_continuous(labels = comma, limits = c(0, 150), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_size_continuous(labels = comma, range = c(0.4, 2.5), guide = "none") +
  labs(
    x = "Distance between sites (kb)",
    y = expression(mean~r^2),
    title = "LD decay from 0 to 150 kb",
    subtitle = "Pairwise r2, summarized in 1 kb distance bins"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p
ggsave(
  file.path(out_dir, "ld_decay_0_1000kb.png"),
  p,
  width = 8.5,
  height = 5.2,
  dpi = 300
)

message("Wrote LD decay summary and plot to: ", normalizePath(out_dir))
