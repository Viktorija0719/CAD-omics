#!/usr/bin/env Rscript

# DESCRIPTION:
# This script reads a count matrix of miRNA expression data, performs TMM normalization,
# applies log2(CPM + 1) transformation, and writes the normalized matrix to an output file.
#
# USAGE:
# module load gnu8/8.3.0
# module load R/R-4.2.0
# Rscript normalize_tmm_log2cpm.R input_file.txt output_file.txt

main <- function(argv = commandArgs(trailingOnly = TRUE)) {

  if (length(argv) != 2) {
    stop(" Usage: Rscript normalize_tmm_log2cpm.R <input_file> <output_file>")
  }

  input_path <- argv[1]
  output_path <- argv[2]

  suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
  })

  # Read input file
  countTable <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Extract miRNA IDs and counts
  miRNA_ids <- countTable[[1]]
  counts <- countTable[, -1]

  # Normalize and transform
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  norm_matrix <- log2(cpm(dge) + 1)

  # Combine with miRNA names
  output_df <- data.frame(mirna = miRNA_ids, norm_matrix, check.names = FALSE)

  # Write output
  write.table(output_df, file = output_path, sep = ",", quote = FALSE, row.names = FALSE)
}

main()