#!/usr/bin/env Rscript

# ####
# DESCRIPTION:
# This script preprocesses lipid multi-omics data by aligning miRNA, SNP, and clinical datasets
# according to the status file, ensuring all samples in the status file are included.
# The clinical dataset is split into continuous and binary subsets, and near-zero variance features
# are removed from SNP and clinical layers. The output is saved as both RData and RDS objects.

# USAGE (bash):
# Rscript lipid_multiomics_preprocessing.R <mirna_csv> <snp_csv> <clinical_csv> <status_csv> <output_rdata> <output_rds>
# ####

# --- Load required library ---
library(caret)  # for nearZeroVar()

# --- Read arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("USAGE: Rscript R_multi_omics.R <mirna_csv> <snp_csv> <clinical_csv> <status_csv> <output_rdata> <output_rds>")
}
mirna_path    <- args[1]
snp_path      <- args[2]
clinical_path <- args[3]
status_path   <- args[4]
output_rdata  <- args[5]
output_rds    <- args[6]

# --- Load datasets ---
mirna    <- read.csv(mirna_path,    row.names = 1, check.names = FALSE)
snp      <- read.csv(snp_path,      row.names = 1, check.names = FALSE)
clinical <- read.csv(clinical_path, row.names = 1, check.names = FALSE)
status   <- read.csv(status_path,   row.names = 1, check.names = FALSE)

# --- Define continuous and binary clinical variables ---
clinical_continuous <- c(
  "PRS", "Age", "Weight", "Height", "BMI", "Pack_years", "WBC", "RBC", "HGB", "HCT", "PLT", 
  "NE_abs", "LY_abs", "MO_abs", "Eo_Abs", "Total_cholesterol", "LDL", "HDL", "TG", 
  "Creatinine", "ALAT", "ASAT", "Glucose", "HbA1c", "LAD_stenosis", "LCX_stenosis", "RCA_stenosis", 
  "Target_fibrotic", "Target_lipidic", "Target_necrotic", "Target_Calcific", "Plaque_necrolipidic_tissue", 
  "Dist_fibrotic", "Dist_lipidic", "Dist_necrotic", "Dist_calcific"
)

clinical_binary <- c(
  "Sex", "Smoking_1", "Smoking_2", "Positive_family_history", 
  "Art_hipert", "Congestive_heart_failure", "Previous_PCI", "Previous_MI"
)

# --- Use all samples from status ---
samples_all <- rownames(status)

# --- Function to align and expand datasets ---
expand_to_samples <- function(df, samples) {
  df <- df[, intersect(samples, colnames(df)), drop = FALSE]
  missing <- setdiff(samples, colnames(df))
  if (length(missing) > 0) {
    for (s in missing) df[, s] <- NA
  }
  df <- df[, samples, drop = FALSE]
  return(df)
}

mirna <- expand_to_samples(mirna, samples_all)
snp   <- expand_to_samples(snp, samples_all)

clinical <- clinical[match(intersect(samples_all, clinical$sample), clinical$sample), , drop = FALSE]
missing_samples_clin <- setdiff(samples_all, clinical$sample)
if (length(missing_samples_clin) > 0) {
  for (s in missing_samples_clin) clinical[nrow(clinical) + 1, ] <- NA
  rownames(clinical)[(nrow(clinical) - length(missing_samples_clin) + 1):nrow(clinical)] <- missing_samples_clin
}
rownames(clinical) <- clinical$sample
clinical$sample <- NULL
clinical <- clinical[samples_all, , drop = FALSE]

# --- Separate continuous and binary clinical layers ---
clinical_cont <- clinical[, intersect(clinical_continuous, colnames(clinical)), drop = FALSE]
clinical_bin  <- clinical[, intersect(clinical_binary, colnames(clinical)), drop = FALSE]

# --- Ensure correct orientation ---
mirna_t         <- t(mirna)
snp_t           <- t(snp)
clinical_cont_t <- clinical_cont
clinical_bin_t  <- clinical_bin

# --- Filter near-zero variance features ---
nzv_snp <- nearZeroVar(snp_t)
if (length(nzv_snp) > 0) snp_t <- snp_t[, -nzv_snp, drop = FALSE]

nzv_clin_cont <- nearZeroVar(clinical_cont_t)
if (length(nzv_clin_cont) > 0) clinical_cont_t <- clinical_cont_t[, -nzv_clin_cont, drop = FALSE]

nzv_clin_bin <- nearZeroVar(clinical_bin_t)
if (length(nzv_clin_bin) > 0) clinical_bin_t <- clinical_bin_t[, -nzv_clin_bin, drop = FALSE]

# --- Create multi-omics list object ---
my_multiomics <- list(
  data = list(
    mirna           = as.matrix(mirna_t),
    snp             = as.matrix(snp_t),
    clinical_cont   = as.matrix(clinical_cont_t),
    clinical_binary = as.matrix(clinical_bin_t),
    status          = as.factor(status$status)
  )
)

# --- Save outputs ---
save(my_multiomics, file = output_rdata)
saveRDS(my_multiomics, file = output_rds)
message("\nMulti-omics object saved to:\n", output_rdata, "\n", output_rds)
