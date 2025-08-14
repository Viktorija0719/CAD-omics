# --------------------------
# 1. Load your multi-omics dataset
# --------------------------
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



# --------------------------
# 2) Fit iClusterBayes
# --------------------------
# K=1 → 2 clusters (binary setting). You can try K=1:2 later.
set.seed(42)
tic("iClusterBayes")
fit_tune <- tune.iClusterBayes(
  dt1 = data_list$miRNA,
  dt2 = data_list$SNP,
  dt3 = data_list$Clinical_cont,
  dt4 = data_list$Clinical_binary,
  type = c("gaussian","gaussian","gaussian","binomial"),
  K = 2,                # 2 clusters for binary status
  n.burnin = 2000,
  n.draw   = 3000,
  prior.gamma = rep(0.2, 4),
  sdev = 0.5,
  cpus = 4
)
toc()

# 3) Pick best model by BIC
bic_vals <- sapply(fit_tune$fit, function(f) f$BIC)
best_id  <- which.min(bic_vals)
final_model <- fit_tune$fit[[best_id]]

# 4) Clusters
clusters <- final_model$clusters  # 1/2 cluster labels (arbitrary IDs)

# --------------------------
# 5) Map clusters → predicted labels (0/1) for scoring
# --------------------------
truth <- factor(as.character(status), levels = c("0","1"))

# Try both mappings (1->"0"/2->"1") vs (1->"1"/2->"0") and keep the best
map_pred <- function(map01) {
  # map01 is c("1"="0","2"="1") or c("1"="1","2"="0")
  pred <- factor(unname(map01[as.character(clusters)]), levels = c("0","1"))
  cm <- table(truth, pred)
  acc <- sum(diag(cm)) / sum(cm)
  list(pred = pred, cm = cm, acc = acc)
}
alt1 <- map_pred(c("1"="0","2"="1"))
alt2 <- map_pred(c("1"="1","2"="0"))
pick <- if (alt1$acc >= alt2$acc) alt1 else alt2
pred_lab <- pick$pred
cm       <- pick$cm

cat("Confusion matrix (mapped clusters):\n"); print(cm)
cat("Accuracy:", round(pick$acc, 3), "\n")

# --------------------------
# 6) Same metrics as DIABLO: Accuracy, MCC, MAE, PRF, macro-F1
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



perf_icb <- PerformanceAnalysisGeneral(truth, pred_lab)
cat("ICB Accuracy:", perf_icb$accuracy,
    "MCC:", perf_icb$mcc,
    "macro-F1:", perf_icb$macro_F1, "\n")

# Label-invariant clustering scores
ARI <- mclust::adjustedRandIndex(truth, clusters)
NMI <- aricode::NMI(truth, clusters)
cat("Label-invariant metrics — ARI:", round(ARI,3), " NMI:", round(NMI,3), "\n")

# --------------------------
# 7) Save everything in one file
# --------------------------



saveRDS(list(
  final_model = final_model,
  perf_icb = perf_icb
), file = "/home/vik/RSU_work/CAD-omics/results/iClusterBayes_CADomics_results.RDS")

