#!/usr/bin/env Rscript

# DESCRIPTION:
#   PCA of miRNA expression with points colored by a 2-class status.
#   miRNA file: rows = miRNAs, columns = samples (col 1 = miRNA/feature ID).
#   Status file: CSV with columns {sample,status}, where status is 0 or 1.
#   Group labels for 0/1 and optional plot output path are supplied on the command line.
#
# USAGE:
#   Rscript PCA_plot.R <mirna_csv> <status_csv> <group0_label> <group1_label> [output_plot]
#
# EXAMPLE:
#   Rscript PCA_plot.R \
#     data/processed/mirna/lipid_tmm_normalized_adj_sex_age_data.csv \
#     data/processed/status/lipid_status.csv \
#     "Non lipid rich" "Lipid rich" \
#     plots/lipid_PCA.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

# ---- Args ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop("USAGE: Rscript PCA_plot.R <mirna_csv> <status_csv> <group0_label> <group1_label> [output_plot]")
}
mirna_path   <- args[1]
status_path  <- args[2]
label_group0 <- args[3]
label_group1 <- args[4]
output_plot  <- ifelse(length(args) == 5, args[5], "")

cat("[INFO] miRNA:", mirna_path, "\n")
cat("[INFO] Status:", status_path, "\n")
cat("[INFO] Labels: 0 ->", label_group0, "| 1 ->", label_group1, "\n")
if (nzchar(output_plot)) cat("[INFO] Output plot:", output_plot, "\n")

# ---- Load miRNA ----
dat <- read.csv(mirna_path, check.names = FALSE)
if (ncol(dat) < 3) stop("miRNA file seems to have too few columns. Expect: 1st = miRNA ID, others = samples.")
expr_mat <- as.matrix(dat[, -1, drop = FALSE])
rownames(expr_mat) <- dat[[1]]
dat2 <- t(expr_mat)

# ---- Load status ----
status_df <- read.csv(status_path, stringsAsFactors = FALSE, check.names = FALSE)
colnames(status_df) <- trimws(tolower(colnames(status_df)))
if (!all(c("sample", "status") %in% colnames(status_df))) {
  stop("Status file must contain columns: 'sample' and 'status'")
}
status_df$status <- factor(status_df$status, levels = c(0, 1), labels = c(label_group0, label_group1))

# ---- Align samples ----
common <- intersect(rownames(dat2), status_df$sample)
if (length(common) < 3) stop("Not enough overlapping samples for PCA.")
dat2 <- dat2[common, , drop = FALSE]
status_df <- status_df[match(common, status_df$sample), , drop = FALSE]

# ---- PCA ----
dat.pca <- prcomp(dat2, center = TRUE, scale. = FALSE)
explained_variance <- summary(dat.pca)$importance[2, 1:2] * 100

# ---- Plot df ----
pca_df <- data.frame(
  Sample = rownames(dat2),
  PC1 = dat.pca$x[, 1],
  PC2 = dat.pca$x[, 2],
  status = status_df$status,
  stringsAsFactors = FALSE
)

# ---- Aesthetics ----
shape_vals <- setNames(c(21, 21), c(label_group0, label_group1))
col_vals   <- setNames(c("seagreen", "firebrick"), c(label_group0, label_group1))
fill_vals  <- setNames(c(alpha("seagreen", 0.5), alpha("firebrick", 0.5)), c(label_group0, label_group1))

# ---- Plot ----
plt <- ggplot(pca_df, aes(x = PC1, y = PC2, color = status, shape = status, fill = status)) +
  geom_point(size = 2, stroke = 1.2) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 100) +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(values = col_vals) +
  scale_fill_manual(values = fill_vals) +
  stat_ellipse(type = "t", linetype = "dashed", linewidth = 0.5) +
  theme_minimal() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  labs(
    title = sprintf("PCA by %s vs %s", label_group0, label_group1),
    x = sprintf("PC1 (%.2f%% variance)", explained_variance[1]),
    y = sprintf("PC2 (%.2f%% variance)", explained_variance[2])
  )

# ---- Output ----
if (nzchar(output_plot)) {
  # ensure directory exists
  out_dir <- dirname(output_plot)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(output_plot, plt, width = 7, height = 5, dpi = 300)
  cat("[INFO] Saved plot to:", output_plot, "\n")
} else {
  print(plt)
}
