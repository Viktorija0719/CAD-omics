#!/usr/bin/env Rscript

# DESCRIPTION:
#   Filter limma differential expression result CSVs in a given directory by
#   log2 fold-change and adjusted p-value thresholds.
#
# INPUT:
#   Directory containing limma_DEGs_*.csv files.
#
# FILTERS (defaults):
#   |log2FC| > 1  AND  adjusted p-value < 0.05
#
# OUTPUT:
#   For each "<name>.csv", writes "<name>_filtered.csv" to the same directory.
#
# USAGE:
#   docker run -it --rm \
#     -v $(pwd):/app \
#     cad-omics \
#     Rscript scripts/filter_limma_results.R results/stability_limma_adj_sex_age_data
#
# DEPENDENCIES: base R only

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript filter_limma_results.R <results_dir>")
results_dir <- args[1]
if (!dir.exists(results_dir)) stop("Directory not found: ", results_dir)

# --- Helpers ---
normalize_names <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}
find_col <- function(cols, patterns) {
  cn <- normalize_names(cols)
  for (p in patterns) {
    hit <- which(grepl(p, cn, perl = TRUE))
    if (length(hit) >= 1) return(cols[hit[1]])
  }
  NA_character_
}

# --- Filtering function ---
filter_limma_results <- function(input_path,
                                 logFC_threshold = log2(2),
                                 pval_threshold  = 0.05) {

  dat <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)

  logfc_col <- find_col(colnames(dat), c("logfc", "log2fc"))
  fdr_col   <- find_col(colnames(dat), c("adjpval", "adjpvalue", "fdr", "padj", "qvalue"))

  if (is.na(logfc_col)) stop("No logFC column found in: ", input_path)
  if (is.na(fdr_col))   stop("No adjusted p-value column found in: ", input_path)

  dat[[logfc_col]] <- suppressWarnings(as.numeric(dat[[logfc_col]]))
  dat[[fdr_col]]   <- suppressWarnings(as.numeric(dat[[fdr_col]]))

  keep <- (abs(dat[[logfc_col]]) > logFC_threshold) & (dat[[fdr_col]] < pval_threshold)
  filtered <- dat[keep, , drop = FALSE]

  out_file <- sub("\\.csv$", "_filtered.csv", basename(input_path), ignore.case = TRUE)
  out_path <- file.path(dirname(input_path), out_file)

  write.csv(filtered, out_path, row.names = FALSE)
  cat("Saved:", out_path, "(", nrow(filtered), "rows)\n")
}

# --- Run over all matching files ---
files <- list.files(results_dir, pattern = "^limma_DEGs_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No limma_DEGs_*.csv files found in ", results_dir)

for (f in files) filter_limma_results(f)
