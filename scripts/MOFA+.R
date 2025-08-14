#############################
## MOFA+ on CAD dataset   ##
#############################

suppressPackageStartupMessages({
  library(MOFA2)
  library(ggplot2)
})


my_multiomics <- readRDS("/home/vik/RSU_work/CAD-omics/data/processed/lipid_multiomics_imputed_binary.rds")

# Extract omics blocks
data_list <- list(
  miRNA           = as.matrix(my_multiomics$data$mirna),
  SNP             = as.matrix(my_multiomics$data$snp),
  Clinical_cont   = as.matrix(my_multiomics$data$clinical_cont),
  Clinical_binary = as.matrix(my_multiomics$data$clinical_binary)
)

# --------------------------
# 2. Clean and match samples
# --------------------------
# Filter out samples with all NA per block
data_list <- lapply(data_list, function(x) x[complete.cases(x), , drop = FALSE])

# Get intersecting sample names
common_samples <- Reduce(intersect, lapply(data_list, rownames))

# Subset data blocks
data_list <- lapply(data_list, function(x) x[common_samples, , drop = FALSE])

# Prepare sample labels (status)
status <- my_multiomics$data$status
names(status) <- rownames(my_multiomics$data$mirna)
status <- status[common_samples]
status <- droplevels(status)
status_labels <- as.character(status)

# Check dimensions
cat("Sample sizes per block:\n")
print(lapply(data_list, dim))
cat("Status table:\n")
print(table(status))




# MOFA expects features x samples. If you currently have samples x features, transpose.
transpose_if_needed <- function(M) {
  if (!is.null(rownames(M)) && nrow(M) == length(common_samples)) {
    # likely samples x features -> transpose to features x samples
    TM <- t(M)
    colnames(TM) <- rownames(M)
    rownames(TM) <- colnames(M)
    return(TM)
  } else {
    return(M)
  }
}
data_list <- lapply(data_list, transpose_if_needed)

# --------------------------
# 3) Build and train MOFA+
# --------------------------
views <- list(
  miRNA               = data_list$miRNA,            # continuous
  SNPs                = data_list$SNP,              # 0/1/2 dosage -> gaussian is fine
  Clinical_continuous = data_list$Clinical_cont,    # continuous
  Clinical_binary     = data_list$Clinical_binary   # binary (0/1)
)

# Create object
MO <- create_mofa(views)

# Model options
model_opts <- get_default_model_options(MO)
model_opts$num_factors <- 10
model_opts$likelihoods <- c(
  miRNA               = "gaussian",
  SNPs                = "bernoulli",   # keep as gaussian for dosage 0/1/2
  Clinical_continuous = "gaussian",
  Clinical_binary     = "bernoulli"
)

# Data & training options
data_opts <- get_default_data_options(MO)
train_opts <- get_default_training_options(MO)
train_opts$seed <- 600
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "medium"

# Prepare & run
MO <- prepare_mofa(
  MO,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
M <- run_mofa(MO, use_basilisk = TRUE)




# # before run_mofa()
# options(timeout = 600)                        # give it 10 min instead of 60s
# options(download.file.method = "libcurl")     # more robust on Linux
# # optional: use wget if you have it
# # options(download.file.method = "wget")
# 
# # put the cache somewhere you can write and reuse across runs
# Sys.setenv(BASILISK_CACHE_DIR = "~/.cache/basilisk")
# Sys.setenv(BASILISK_DOWNLOAD_DIR = "~/.cache/basilisk")
# 
# # (if you are behind a proxy, set these)
# # Sys.setenv(https_proxy="http://user:pass@host:port", http_proxy="http://user:pass@host:port")
# 
# # also, pick an explicit output file so it doesn’t go to /tmp
# M <- run_mofa(MO, use_basilisk = TRUE,
#               outfile = "/home/vik/RSU_work/CAD-omics/results/mofa_model.hdf5")





# Attach metadata for convenience
samples_metadata(M) <- data.frame(
  sample = common_samples,
  status = status,
  row.names = common_samples
)

# Factors matrix (samples x factors) from the (single) group
F_mat <- get_factors(M)[[1]]  # samples x k

# --------------------------
# 4) Stratified 75/25 split
# --------------------------
set.seed(600)
p <- 0.75
train_idx <- integer(0)
for (lev in levels(status)) {
  idx <- which(status == lev)
  n_take <- max(1, floor(p * length(idx)))
  train_idx <- c(train_idx, sample(idx, n_take))
}
train_idx <- sort(unique(train_idx))
test_idx  <- setdiff(seq_along(status), train_idx)

Y_train <- droplevels(status[train_idx])
Y_test  <- droplevels(status[test_idx])

F_train <- F_mat[train_idx, , drop = FALSE]
F_test  <- F_mat[test_idx, , drop = FALSE]

# Safety: ensure both classes appear in train
if (length(unique(Y_train)) < 2L) stop("Training set lost a class. Table: ", paste(names(table(Y_train)), table(Y_train), collapse="; "))

# --------------------------
# 5) Classifier on factors
# --------------------------




## 5) Classifier on factors (logistic regression)
# Ensure factor matrices share identical column names
stopifnot(identical(colnames(F_train), colnames(F_test)))

df_train <- data.frame(y = as.numeric(Y_train) - 1, F_train)
fit <- stats::glm(y ~ ., data = df_train, family = stats::binomial())

# Predict on test (use stats::predict explicitly)
df_test <- as.data.frame(F_test)
probs   <- stats::predict(fit, newdata = df_test, type = "response")
preds   <- ifelse(probs >= 0.5, "1", "0")
preds   <- factor(preds, levels = levels(status))  # keep same label levels

# Quick sanity check
table(Pred = preds, Truth = Y_test)

# 
# # Logistic regression (simple & interpretable)
# df_train <- data.frame(y = as.numeric(Y_train) - 1, F_train)
# fit <- glm(y ~ ., data = df_train, family = binomial())
# 
# # Predict on test
# probs <- predict(fit, newdata = data.frame(F_test), type = "response")
# 
# 
# 
# probs <- predict(fit, newdata = data.frame(F_test))
# preds <- ifelse(probs >= 0.5, "1", "0")
# preds <- factor(preds, levels = levels(status))

# --------------------------
# 6) Metrics: Confusion, Accuracy, MCC, MAE, per-class PRF, macro-F1
# --------------------------









PerformanceAnalysisGeneral <- function(truth, predicted) {
  truth <- factor(as.character(truth), levels = c("0", "1"))
  predicted <- factor(as.character(predicted), levels = c("0", "1"))
  cm <- table(truth, predicted)
  
  # Accuracy
  acc <- sum(diag(cm)) / sum(cm)
  
  # MCC (binary)
  if (all(dim(cm) == c(2, 2))) {
    TP <- cm[2, 2]
    TN <- cm[1, 1]
    FP <- cm[1, 2]
    FN <- cm[2, 1]
    denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
    mcc <- if (denom > 0) (TP * TN - FP * FN) / sqrt(denom) else NA_real_
  } else {
    mcc <- NA_real_
  }
  
  # Precision, Recall, F1 per class
  precision <- diag(cm) / colSums(cm)
  recall    <- diag(cm) / rowSums(cm)
  f1        <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  
  per_class <- data.frame(
    class     = levels(truth),
    precision = round(precision, 3),
    recall    = round(recall, 3),
    F1        = round(f1, 3)
  )
  
  list(
    confusion = cm,
    accuracy  = round(acc, 3),
    mcc       = round(mcc, 3),
    per_class = per_class,
    macro_F1  = round(mean(f1), 3)
  )
}

perf_mofa <- PerformanceAnalysisGeneral(Y_test, preds)

cat("\n=== MOFA+ classifier (logit on factors) — test set ===\n")
print(perf_mofa$confusion)
cat(
  "Accuracy:", perf_mofa$accuracy,
  "| MCC:", perf_mofa$mcc,
  "| macro-F1:", perf_mofa$macro_F1, "\n"
)

# --------------------------
# 7) Save everything in one file
# --------------------------
save_path <- "/home/vik/RSU_work/CAD-omics/results/mofa2_cad_omics_allinone.RDS"
saveRDS(list(
  MOFA_object   = MO,
  MOFA_model    = M,
  factors       = F_mat,          # samples x factors
  train_index   = train_idx,
  test_index    = test_idx,
  status        = status,
  glm_fit       = fit,
  probs_test    = probs,
  preds_test    = preds,
  metrics       = perf_mofa
), file = save_path)
cat("Saved results to:", save_path, "\n")

# (Optional) A few QC plots you can call interactively:
# plot_variance_explained(M, plot_total = TRUE)
# plot_factor(M, factors = 1, color_by = "status", add_violin = TRUE, dodge = TRUE)
# plot_factor_cor(M)


