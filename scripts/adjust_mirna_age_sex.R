#!/usr/bin/env Rscript

# --- DESCRIPTION ---
# Adjusts miRNA expression data for Age and Sex using linear regression.
# 
# INPUTS (specified as arguments when calling from bash):
#   1. miRNA CSV file (rows = miRNAs, columns = samples)
#   2. Metadata CSV file (must contain columns: sample, Age, Sex)
#   3. Output CSV path (adjusted miRNA values will be saved here)
#
# USAGE EXAMPLE:
#   Rscript adjust_mirna_age_sex.R <mirna_csv> <meta_csv> <output_csv>
#
# EXAMPLE:
#   Rscript adjust_mirna_age_sex.R \
#       data/processed/mirna/stability_tmm_normalized_data.csv \
#       data/processed/clinical/stability_meta_missing.csv \
#       data/processed/mirna/stability_tmm_normalized_adj_sex_age_data.csv

# --- Read arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("USAGE: Rscript adjust_mirna_age_sex.R <mirna_csv> <meta_csv> <output_csv>")
}
mirna_path  <- args[1]
meta_path   <- args[2]
output_path <- args[3]

cat("\n[INFO] miRNA file:", mirna_path)
cat("\n[INFO] Metadata file:", meta_path)
cat("\n[INFO] Output file:", output_path, "\n")

# --- Load data ---
mirna <- read.csv(mirna_path, row.names = 1, check.names = FALSE)
meta  <- read.csv(meta_path, check.names = FALSE)

# --- Remove any empty column names ---
meta <- meta[, !(names(meta) == "" | is.na(names(meta)))]

# --- Clean names ---
names(meta) <- trimws(names(meta))

# --- Drop index column if present ---
if (!"sample" %in% names(meta) & names(meta)[1] == "") {
  names(meta)[1] <- "row_index"
}
if ("row_index" %in% names(meta)) {
  meta$row_index <- NULL
}

# --- Extract sample names ---
mirna_samples <- colnames(mirna)
meta_samples  <- meta$sample

# --- Keep intersecting samples ---
common_samples <- intersect(mirna_samples, meta_samples)

# --- Align order ---
mirna <- mirna[, common_samples]
meta  <- meta[match(common_samples, meta$sample), ]

# --- Check alignment ---
stopifnot(rownames(t(mirna)) == meta$sample)

# --- Transpose for modeling ---
mirna_t <- t(mirna)

# --- Ensure correct types ---
meta$Sex <- as.factor(meta$Sex)

# --- Adjust for Age and Sex ---
adjusted_mirna <- apply(mirna_t, 2, function(y) {
  residuals(lm(y ~ Age + Sex, data = meta))
})

# --- Back to miRNA-style format ---
adjusted_mirna <- t(adjusted_mirna)
rownames(adjusted_mirna) <- rownames(mirna)
colnames(adjusted_mirna) <- colnames(mirna)

# --- Save ---
write.csv(adjusted_mirna, output_path)
cat("\n✅ Adjusted miRNA data saved to:", output_path, "\n")
