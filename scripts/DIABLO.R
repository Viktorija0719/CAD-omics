if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mixOmics", force = TRUE)


rm(list=ls())
######################
### LOAD LIBRARIES ###
######################

library(mixOmics)
library(tictoc)
library(ggplot2)
library(dplyr)
library(reshape2)
# Optional parallel (for tune/perf): BiocParallel::register(BiocParallel::MulticoreParam(2))

# set.seed(600)
# 
# ##########################
# ### LOAD YOUR DATASET  ###
# ##########################
# 
# # Expect: my_multiomics$data list contains:
# #   mirna, snp, clinical_cont, clinical_binary, status
# my_multiomics <- readRDS("/home/vik/RSU_work/CAD-omics/data/processed/stability_multiomics_imputed_binary.rds")
# 
# 
# # Build blocks
# clinical_combined <- cbind(my_multiomics$data$clinical_cont,
#                            my_multiomics$data$clinical_binary)
# 
# data_blocks <- list(
#   miRNA    = my_multiomics$data$mirna,
#   SNP      = my_multiomics$data$snp,
#   Clinical = clinical_combined
# )
# 
# # Outcome (binary factor with levels "0","1")
# Y <- as.factor(my_multiomics$data$status)
# 
# # Ensure all blocks are numeric matrices
# to_numeric <- function(M){ M <- as.matrix(M); storage.mode(M) <- "double"; M }
# data_blocks <- lapply(data_blocks, to_numeric)
# 
# # Align sample IDs across blocks (by intersected rownames)
# common_ids <- Reduce(intersect, lapply(data_blocks, rownames))
# stopifnot(length(common_ids) > 0)
# data_blocks <- lapply(data_blocks, function(M) M[common_ids, , drop = FALSE])
# 
# # Align labels to blocks using sample IDs if names exist; otherwise assume same order
# if (!is.null(names(Y))) {
#   Y <- Y[common_ids]
# } else {
#   if (length(Y) != length(common_ids)) stop("Y length does not match sample count.")
#   names(Y) <- rownames(data_blocks$miRNA)
#   Y <- Y[common_ids]
# }
# Y <- droplevels(as.factor(as.character(Y)))
# stopifnot(!any(is.na(Y)))
# stopifnot(all(levels(Y) %in% c("0","1")))
# 
# cat("Block dims:\n"); print(lapply(data_blocks, dim))
# cat("Class table:\n"); print(table(Y))








set.seed(600)


my_multiomics <- readRDS("/home/vik/RSU_work/CAD-omics/data/processed/stability_multiomics_imputed_binary.rds")

# Build blocks
clinical_combined <- cbind(my_multiomics$data$clinical_cont,
                           my_multiomics$data$clinical_binary)

data_blocks <- list(
  miRNA    = my_multiomics$data$mirna,
  SNP      = my_multiomics$data$snp,
  Clinical = clinical_combined
)



# --------------------------
# 2. Clean and match samples
# --------------------------
# Filter out samples with all NA per block
data_blocks <- lapply(data_blocks, function(x) x[complete.cases(x), , drop = FALSE])

# Get intersecting sample names
common_samples <- Reduce(intersect, lapply(data_blocks, rownames))

# Subset data blocks
data_blocks <- lapply(data_blocks, function(x) x[common_samples, , drop = FALSE])

# Prepare sample labels (status)
status <- my_multiomics$data$status
names(status) <- rownames(my_multiomics$data$mirna)
status <- status[common_samples]
status <- droplevels(status)
status_labels <- as.character(status)

# Check dimensions
cat("Sample sizes per block:\n")
print(lapply(data_blocks, dim))
cat("Status table:\n")
print(table(status))




# Outcome labels aligned AFTER filtering
status <- my_multiomics$data$status
names(status) <- rownames(my_multiomics$data$mirna)
status <- status[common_samples]
Y <- droplevels(as.factor(as.character(status)))  # final Y






####################################
### GLOBAL SETTINGS / SMALL UTILS ###
####################################

maxComp  <- 3   # upper bound for selecting number of components
dist_metric <- "centroids.dist"  # more stable than mahalanobis on small N
# keepX grid helper (clipped by block size)
keep_grid <- function(M){
  base <- c(5:9, seq(10, 18, 2), seq(20, 30, 5))
  unique(pmin(base, ncol(M)))
}

# Drop exact zero-variance columns
rm_zerovar <- function(M){
  v <- apply(M, 2, function(x) var(x, na.rm = TRUE))
  keep <- which(!is.na(v) & v > 0)
  if(length(keep) == 0) stop("All columns removed due to zero variance in a block.")
  M[, keep, drop = FALSE]
}

################################
### PERFORMANCE HELPER ###
################################



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

#############################
### CORE ANALYSIS FUNCTION ###
#############################

mixAnalysis <- function(blocks, class_labels, dataname){
  # 1) Stratified 75/25 split
  class_labels <- factor(class_labels)
  tabY <- table(class_labels)
  if (any(tabY < 2)) stop("Each class needs at least 2 samples: ", paste(names(tabY), tabY, collapse=" "))
  
  set.seed(600)
  p <- 0.75
  train_ind <- integer(0)
  for (lev in levels(class_labels)) {
    idx <- which(class_labels == lev)
    n_take <- max(1, floor(p * length(idx)))
    train_ind <- c(train_ind, sample(idx, n_take))
  }
  train_ind <- sort(unique(train_ind))
  
  # Split
  data_train <- lapply(blocks, function(M) M[train_ind, , drop = FALSE])
  data_test  <- lapply(blocks, function(M) M[-train_ind, , drop = FALSE])
  Y_train <- droplevels(class_labels[train_ind])
  Y_test  <- droplevels(class_labels[-train_ind])
  
  if (nlevels(Y_train) < 2) stop("Training set lost a class.")
  if (nlevels(Y_test)  < 2) message("Warning: test set lost a class (small sample).")
  
  # 2) Filter zero-variance on TRAIN only; align TEST columns
  data_train <- lapply(data_train, rm_zerovar)
  data_test  <- Map(function(te, tr) te[, colnames(tr), drop = FALSE], data_test, data_train)
  
  # 3) Scale by TRAIN stats
  centers <- list(); scales <- list()
  data_train <- lapply(seq_along(data_train), function(i){
    Z <- scale(data_train[[i]])
    centers[[i]] <<- attr(Z, 'scaled:center')
    scales[[i]]  <<- attr(Z, 'scaled:scale')
    Z
  })
  names(data_train) <- names(blocks)
  data_test <- lapply(seq_along(data_test), function(i){
    scale(data_test[[i]], center = centers[[i]], scale = scales[[i]])
  })
  names(data_test) <- names(blocks)
  
  tictoc::tic(dataname)
  
  # 4) Design matrix
  design <- matrix(0.1, ncol = length(data_train), nrow = length(data_train),
                   dimnames = list(names(data_train), names(data_train)))
  diag(design) <- 0
  
  # 5) Choose ncomp via quick CV (cap folds by smallest class)
  folds_used <- min(5, as.integer(min(table(Y_train))))
  if (folds_used < 2) folds_used <- 2
  
  base_model <- block.splsda(X = data_train, Y = Y_train, ncomp = maxComp,
                             design = design, near.zero.var = TRUE)
  perf_base  <- perf(base_model, validation = 'Mfold', folds = folds_used,
                     nrepeat = 3, dist = dist_metric)
  ER <- "Overall.ER"
  #ncomp <- perf_base$choice.ncomp$WeightedVote[ER, dist_metric]
  #if (is.na(ncomp) || ncomp < 1) ncomp <- 1
  ncomp <- 2
  
  # 6) Tune keepX per block
  test.keepX <- lapply(data_train, keep_grid)
  tune_res <- tune.block.splsda(X = data_train, Y = Y_train, ncomp = ncomp,
                                test.keepX = test.keepX, design = design,
                                validation = 'Mfold', folds = folds_used,
                                nrepeat = 1, dist = dist_metric)
  list.keepX <- lapply(tune_res$choice.keepX, function(x) x[ncomp])
  
  # 7) Final model
  final_model <- block.splsda(X = data_train, Y = Y_train, ncomp = ncomp,
                              keepX = list.keepX, design = design,
                              near.zero.var = TRUE)
  
  # 8) Training CV summary (optional)
  perf_final <- perf(final_model, validation = 'Mfold', folds = folds_used,
                     nrepeat = 3, dist = dist_metric)
  
  # 9) Predict held-out test set
  pred <- predict(final_model, newdata = data_test)
  #pred <- predict(final_model, data_test)
  #pred <- mixOmics::predict(final_model, newdata = data_test)
  
  preds_class <- pred$WeightedVote[[dist_metric]][, ncomp, drop = TRUE]
  
  tictoc::toc()
  
  list(
    Y_test      = Y_test,
    pred        = pred,
    preds_class = preds_class,
    ncomp       = ncomp,
    list.keepX  = list.keepX,
    perf        = perf_final,
    final_model = final_model
  )
}








##########################
### RUN WITH YOUR DATA ###
##########################


res <- mixAnalysis(data_blocks, Y, dataname = "CAD_omics")

# Confusion matrix + metrics
cm <- get.confusion_matrix(truth = res$Y_test, predicted = res$preds_class)
print(cm)
acc <- sum(diag(cm)) / sum(cm)
cat("Test Accuracy:", round(acc, 3), "\n")

perf_test <- PerformanceAnalysisGeneral(res$Y_test, res$preds_class)
perf_test

##############################
### OPTIONAL: KEY PLOTS    ###
##############################
# final <- res$final_model; ncomp <- res$ncomp
# plotDiablo(final, ncomp = min(2, ncomp), legend = TRUE, title = "DIABLO global map")
# plotIndiv(final, ncomp = min(2, ncomp), ind.names = FALSE, legend = TRUE, title = "Samples per block")
# plotArrow(final, ncomp = min(2, ncomp), ind.names = FALSE, legend = TRUE, title = "Arrow plot")
# plotVar(final, comp = 1:min(2, ncomp), var.names = FALSE, style = 'graphics', legend = TRUE,
#         pch = c(16,17,15), cex = c(1.8,1.8,1.8),
#         col = c('darkorchid','steelblue','seagreen'))
# circosPlot(final, cutoff = 0.6,
#            color.blocks = c('darkorchid','steelblue','seagreen'),
#            color.cor = c("chocolate3","grey20"))
# network(final, blocks = c(1,2,3),
#         color.node = c('darkorchid','steelblue','seagreen'), cutoff = 0.4)
# plotLoadings(final, block = "miRNA",   comp = 1, contrib = "max", method = "median")
# plotLoadings(final, block = "SNP",     comp = 1, contrib = "max", method = "median")
# plotLoadings(final, block = "Clinical",comp = 1, contrib = "max", method = "median")
# cimDiablo(final)


saveRDS(list(
  res       = res,
  perf_test = perf_test,
  confusion = cm
), file = "/home/vik/RSU_work/CAD-omics/results/stability_DIABLO_common_binary_imputed.RDS")



