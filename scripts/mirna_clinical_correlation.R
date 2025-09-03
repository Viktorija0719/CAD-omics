#!/usr/bin/env Rscript

# DESCRIPTION:
#   Correlate normalized miRNA expression with clinical variables.
#   Inputs:
#     (1) Normalized miRNA matrix CSV: rows = miRNAs, col1 = miRNA ID, other columns = samples.
#     (2) Clinical CSV: rows = samples, must contain 'sample' column and clinical variables.
#   The script aligns by sample ID and analyzes only samples present in the miRNA data.
#
# ANALYSIS:
#   - For continuous clinical variables: Kendall's tau vs each miRNA.
#   - For binary categorical variables: point-biserial correlation vs each miRNA (by coding 0/1).
#   - Benjamini–Hochberg (BH) correction for multiple testing per method.
#
# OUTPUT:
#   One CSV with columns: Variable1 (clinical), Variable2 (miRNA), Correlation Coefficient,
#   p.value, p.adjusted, Method. Filename pattern: <prefix>_combined_results.csv
#
# USAGE (Docker-style):
#   docker run -it --rm \
#     -v $(pwd):/app \
#     cad-omics \
#     Rscript scripts/mirna_clinical_correlation.R \
#       --mirna data/processed/mirna/lipid_tmm_normalized_adj_sex_age_data.csv \
#       --clinical data/processed/clinical/stability_meta_missing.csv \
#       --output_dir results/correlation \
#       --prefix stability_miRNA_corr
#
# DEPENDENCIES: dplyr, readr, stringr, tidyr, optparse

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(optparse)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("--mirna"), type = "character", help = "Normalized miRNA CSV (rows=miRNAs, col1=ID, other cols=samples)"),
  make_option(c("--clinical"), type = "character", help = "Clinical CSV with 'sample' column"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./", help = "Output directory [default: %default]"),
  make_option(c("-p", "--prefix"), type = "character", default = "correlation", help = "Prefix for output files [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$mirna) || is.null(opt$clinical)) {
  stop("Please provide --mirna and --clinical paths. Use --help for details.")
}
if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Clinical variable sets (from your specification) ----------
continuous_vars_full <- c(
  "PRS","Age","Weight","Height","BMI","Pack_years","WBC","RBC","HGB","HCT","PLT",
  "NE_abs","LY_abs","MO_abs","Eo_Abs","Total_cholesterol","LDL","HDL","TG",
  "Creatinine","ALAT","ASAT","Glucose","HbA1c","LAD_stenosis",
  "LCX_stenosis","RCA_stenosis","Target_fibrotic","Target_lipidic","Target_necrotic",
  "Target_Calcific","Plaque_necrolipidic_tissue","Dist_fibrotic","Dist_lipidic","Dist_necrotic","Dist_calcific"
)
categorical_vars_full <- c(
  "Sex","Smoking_1","Smoking_2","Positive_family_history",
  "Art_hipert","Congestive_heart_failure","Previous_PCI","Previous_MI"
)

# ---------- Load data ----------
# miRNA: rows = miRNAs, first col = ID, other cols = samples
mirna_raw <- readr::read_csv(opt$mirna, show_col_types = FALSE)
if (ncol(mirna_raw) < 3) stop("miRNA file must have 1 feature column + >=2 sample columns.")
mirna_ids <- mirna_raw[[1]]
expr_mat  <- as.matrix(mirna_raw[, -1, drop = FALSE])
rownames(expr_mat) <- mirna_ids
suppressWarnings(storage.mode(expr_mat) <- "numeric")
# transpose to samples x miRNAs and add 'sample' column
expr_df <- as.data.frame(t(expr_mat), check.names = FALSE)
expr_df$sample <- rownames(expr_df)

# clinical: must contain 'sample'
clinical <- readr::read_csv(opt$clinical, show_col_types = FALSE)
if (!"sample" %in% names(clinical)) stop("Clinical file must contain a 'sample' column.")

# ---------- Align by samples (keep only those present in miRNA data) ----------
common <- intersect(clinical$sample, expr_df$sample)
if (length(common) < 3) stop("Too few overlapping samples between clinical and miRNA data (need >= 3).")
clinical_f <- clinical %>% filter(sample %in% common)
expr_f     <- expr_df %>% filter(sample %in% common)
merged_data <- clinical_f %>% inner_join(expr_f, by = "sample")

# ---------- Identify miRNA columns ----------
# Accept "hsa-" prefixed (e.g., hsa-miR-...), but if not present, assume all non-clinical columns are miRNAs.
mirna_cols <- names(merged_data)[grepl("^hsa-", names(merged_data))]
if (length(mirna_cols) == 0) {
  # fallback: any columns that were not in clinical plus not 'sample'
  mirna_cols <- setdiff(names(expr_f), "sample")
}
if (length(mirna_cols) == 0) stop("No miRNA columns detected after merge.")

# ---------- Prepare variable lists that are actually present ----------
continuous_vars <- intersect(continuous_vars_full, names(merged_data))
categorical_vars <- intersect(categorical_vars_full, names(merged_data))

# ---------- Helpers ----------
# Kendall's tau for continuous clinical vs miRNA
calculate_kendall_correlation_hsa <- function(data, cont_vars, hsa_vars) {
  if (length(cont_vars) == 0) return(dplyr::tibble())
  # Ensure numerics
  data <- data %>% mutate(dplyr::across(all_of(cont_vars), as.numeric))
  # Expand grid and safe cor.test
  combos <- tidyr::expand_grid(Variable1 = cont_vars, Variable2 = hsa_vars)
  results <- purrr::pmap_dfr(
    combos,
    function(Variable1, Variable2) {
      x <- data[[Variable1]]
      y <- data[[Variable2]]
      # Remove pairs with all NA or <3 non-NA pairs
      ok <- stats::complete.cases(x, y)
      if (sum(ok) < 3) {
        return(tibble::tibble(Variable1 = Variable1, Variable2 = Variable2,
                              Tau = NA_real_, p.value = NA_real_))
      }
      ct <- try(stats::cor.test(x[ok], y[ok], method = "kendall"), silent = TRUE)
      if (inherits(ct, "try-error")) {
        tibble::tibble(Variable1 = Variable1, Variable2 = Variable2,
                       Tau = NA_real_, p.value = NA_real_)
      } else {
        tibble::tibble(Variable1 = Variable1, Variable2 = Variable2,
                       Tau = unname(ct$estimate), p.value = ct$p.value)
      }
    }
  )
  if (nrow(results) == 0) return(results)
  results %>% mutate(p.adjusted = p.adjust(p.value, method = "BH"))
}

# Point-biserial for binary categorical vs miRNA (Pearson after coding 0/1)
calculate_point_biserial_correlation_hsa <- function(data, hsa_vars, disc_vars) {
  out <- list()
  for (disc in disc_vars) {
    # keep only truly binary columns after dropping NA
    vals <- unique(stats::na.omit(data[[disc]]))
    if (length(vals) != 2) next
    # code to 0/1 in a stable order
    binary <- as.integer(factor(data[[disc]], levels = sort(vals))) - 1L
    for (hsa in hsa_vars) {
      x <- binary
      y <- data[[hsa]]
      ok <- stats::complete.cases(x, y)
      if (sum(ok) < 3) next
      ct <- try(stats::cor.test(x[ok], y[ok], method = "pearson"), silent = TRUE)
      if (inherits(ct, "try-error")) next
      out[[paste(disc, hsa, sep = "_")]] <- tibble::tibble(
        Variable1 = disc, Variable2 = hsa,
        Rpb = unname(ct$estimate), p.value = ct$p.value
      )
    }
  }
  df <- dplyr::bind_rows(out)
  if (nrow(df) == 0) return(df)
  df$p.adjusted <- p.adjust(df$p.value, method = "BH")
  df
}

combine_results <- function(kendall, biserial) {
  dplyr::bind_rows(
    kendall %>%
      dplyr::rename(`Correlation Coefficient` = Tau) %>%
      dplyr::mutate(Method = "Kendall"),
    biserial %>%
      dplyr::rename(`Correlation Coefficient` = Rpb) %>%
      dplyr::mutate(Method = "Point-Biserial")
  )
}

# ---------- Run analysis ----------
cat(" Aligning samples and running correlations...\n")
kendall <- calculate_kendall_correlation_hsa(merged_data, continuous_vars, mirna_cols)
biserial <- calculate_point_biserial_correlation_hsa(merged_data, mirna_cols, categorical_vars)
combined <- combine_results(kendall, biserial)

# ---------- Save ----------
out_file <- file.path(opt$output_dir, paste0(opt$prefix, "_combined_results.csv"))
readr::write_csv(combined, out_file)
cat(" Correlation results saved to:", out_file, "\n")
