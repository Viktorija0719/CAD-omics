#!/usr/bin/env Rscript

# DESCRIPTION:
#   Differential expression analysis of normalized miRNA data using limma.
#   Inputs:
#     (1) Normalized miRNA matrix (rows = miRNAs, first column = miRNA ID; columns = samples).
#     (2) Status file (CSV with columns {sample, status}, where 'status' can be 0/1 or >2 groups).
#     (3) Output directory.
#
#   If the status file contains only 0/1 → performs one comparison (Group1 vs Group0).
#   If more than two unique statuses are present → performs all pairwise comparisons
#   (both directions) and outputs one CSV file per contrast.
#
# OUTPUT:
#   For each contrast, a CSV with limma results (logFC, P.Value, adj.P.Val) and
#   median expression for the compared groups.
#
# USAGE:
#   docker run -it --rm \
#     -v $(pwd):/app \
#     cad-omics \
#     Rscript scripts/limma_DE_analysis.R <mirna_file> <status_file.csv> <output_directory>
#
# EXAMPLE:
#   docker run -it --rm \
#     -v $(pwd):/app \
#     cad-omics \
#     Rscript scripts/limma_DE_analysis.R \
#       data/processed/mirna/lipid_tmm_normalized_adj_sex_age_data.csv \
#       data/processed/status/lipid_status.csv \
#       results/limma_output
#
# DEPENDENCIES: limma

suppressPackageStartupMessages({
  library(limma)
})

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) != 3) {
    stop("Usage: Rscript limma_DE_analysis.R <mirna_file> <status_file.csv> <output_directory>")
  }

  mirna_file <- argv[1]
  status_file <- argv[2]
  output_dir  <- argv[3]
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # ---- Load normalized miRNA matrix (first col = ID, the rest = samples) ----
  read_mirna <- function(path) {
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      m <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      m <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    }
    if (ncol(m) < 3) stop("miRNA file must have 1 feature column + >=2 sample columns.")
    expr <- as.matrix(m[, -1, drop = FALSE])
    rownames(expr) <- m[[1]]
    suppressWarnings(storage.mode(expr) <- "numeric")
    expr
  }

  # ---- Load status file (sample, status) ----
  status_df <- read.csv(status_file, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(status_df) <- tolower(trimws(colnames(status_df)))
  if (!all(c("sample","status") %in% colnames(status_df))) {
    stop("Status file must contain columns: 'sample' and 'status'.")
  }

  # ---- Align samples ----
  expr_mat <- read_mirna(mirna_file)
  common <- intersect(colnames(expr_mat), status_df$sample)
  if (length(common) < 3) stop("Not enough overlapping samples between miRNA data and status (need >= 3).")
  expr_mat  <- expr_mat[, common, drop = FALSE]
  status_df <- status_df[match(common, status_df$sample), , drop = FALSE]

  # ---- Build groups from status ----
  levs <- sort(unique(status_df$status))
  group_names <- paste0("Group", seq_along(levs) - 1L)
  status_to_group <- setNames(group_names, levs)
  groups <- factor(status_to_group[as.character(status_df$status)], levels = group_names)

  # ---- Design & fit ----
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(groups)

  logCPM <- expr_mat
  fit <- lmFit(logCPM, design)
  fit$genes <- rownames(logCPM)

  # ---- Build contrasts (fix: do NOT use ':=') ----
  make_pair_name <- function(a, b) sprintf("%svs%s", a, b)

  if (nlevels(groups) == 2) {
    g <- levels(groups)
    contr_vec <- paste0(g[2], "-", g[1])            # e.g., "Group1-Group0"
    contrast_matrix <- makeContrasts(contrasts = contr_vec, levels = design)
    colnames(contrast_matrix) <- make_pair_name(g[2], g[1])  # "Group1vsGroup0"
    comparison_names <- colnames(contrast_matrix)
  } else {
    pairs <- t(combn(levels(groups), 2))
    contr_vec <- character(0)
    comp_names <- character(0)
    for (i in seq_len(nrow(pairs))) {
      A <- pairs[i, 1]; B <- pairs[i, 2]
      contr_vec  <- c(contr_vec,  paste0(A, "-", B), paste0(B, "-", A))
      comp_names <- c(comp_names, make_pair_name(A, B), make_pair_name(B, A))
    }
    contrast_matrix <- makeContrasts(contrasts = contr_vec, levels = design)
    colnames(contrast_matrix) <- comp_names
    comparison_names <- comp_names
  }

  fit2 <- contrasts.fit(fit, contrast_matrix)
  eBayesFit <- eBayes(fit2, trend = TRUE)

  # ---- Medians per group for annotation ----
  group_medians <- lapply(levels(groups), function(g) {
    idx <- which(groups == g)
    apply(logCPM, 1, function(v) median(v[idx], na.rm = TRUE))
  })
  names(group_medians) <- levels(groups)

  add_medians <- function(res, comp_name) {
    parts <- strsplit(comp_name, "vs", fixed = TRUE)[[1]]
    gA <- parts[1]; gB <- parts[2]
    if (gA %in% names(group_medians)) {
      res[[paste0("Median_", gA)]] <- group_medians[[gA]][rownames(res)]
    }
    if (gB %in% names(group_medians)) {
      res[[paste0("Median_", gB)]] <- group_medians[[gB]][rownames(res)]
    }
    res
  }

  # ---- Output per contrast ----
  for (comp in comparison_names) {
    res <- topTable(eBayesFit, coef = comp, number = Inf, sort.by = "none")
    # Add ID as first column (use rownames from topTable)
    res$ID <- rownames(res)
    res <- res[, c("ID", setdiff(names(res), "ID"))]
    res <- add_medians(res, comp)

    out_path <- file.path(output_dir, paste0("limma_DEGs_", comp, ".csv"))
    write.table(res, out_path, sep = ",", quote = FALSE, row.names = FALSE)
    message(" Saved: ", out_path)
  }
}

main()
