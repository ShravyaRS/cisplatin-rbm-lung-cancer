# ================================================================
# Interpretable Rule-Based Prediction of Cisplatin Response
# in Lung Cancer Cell Lines using R.ROSETTA
# ================================================================
# Author: Shravya R S
# Date:  February 16, 2026 
# Data: GDSC2 (Genomics of Drug Sensitivity in Cancer)
# ================================================================

# --- Install required packages (run once) ---
# install.packages(c("BiocManager", "glmnet", "Boruta", "igraph", "devtools"))
# BiocManager::install("PharmacoGx")
# devtools::install_github("komorowskilab/R.ROSETTA")

# --- Load libraries ---
library(PharmacoGx)
library(glmnet)
library(Boruta)
library(R.ROSETTA)
library(igraph)

# ================================================================
# STEP 1: Download GDSC2 dataset
# ================================================================
GDSC2 <- downloadPSet("GDSC_2020(v2-8.2)")

# ================================================================
# STEP 2: Extract cisplatin response for lung cancer
# ================================================================
cell_info <- cellInfo(GDSC2)
lung_cells <- cell_info[cell_info$tissueid == "Lung", ]

sens_data <- summarizeSensitivityProfiles(GDSC2, sensitivity.measure = "aac_recomputed")
cisplatin_response <- sens_data["Cisplatin", ]
cisplatin_response <- cisplatin_response[!is.na(cisplatin_response)]

lung_ids <- rownames(lung_cells)
lung_cisplatin <- cisplatin_response[names(cisplatin_response) %in% lung_ids]

# Classify: median split
threshold <- median(lung_cisplatin)
response_class <- ifelse(lung_cisplatin >= threshold, "Sensitive", "Resistant")

# ================================================================
# STEP 3: Extract gene expression
# ================================================================
rna_data <- summarizeMolecularProfiles(GDSC2, mDataType = "rna",
                                        cell.lines = names(lung_cisplatin))
expr_matrix <- assay(rna_data)
common_cells <- intersect(colnames(expr_matrix), names(lung_cisplatin))
expr_matrix <- expr_matrix[, common_cells]
response_final <- response_class[common_cells]

# ================================================================
# STEP 4: Preprocessing
# ================================================================
# Impute NAs with gene median
for (i in 1:nrow(expr_matrix)) {
  na_idx <- is.na(expr_matrix[i, ])
  if (any(na_idx)) expr_matrix[i, na_idx] <- median(expr_matrix[i, ], na.rm = TRUE)
}

# Remove low-variance genes
gene_vars <- apply(expr_matrix, 1, var)
expr_matrix <- expr_matrix[gene_vars > quantile(gene_vars, 0.25), ]

# Transpose and z-score normalize
expr_norm <- scale(t(expr_matrix))

# ================================================================
# STEP 5: Feature selection - LASSO
# ================================================================
x <- expr_norm
y <- as.factor(response_final)

set.seed(42)
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
coefs <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(coefs)[which(coefs != 0)]
lasso_genes <- lasso_genes[lasso_genes != "(Intercept)"]

# ================================================================
# STEP 6: Feature selection - Boruta
# ================================================================
gene_vars_norm <- apply(expr_norm, 2, var)
top1000 <- names(sort(gene_vars_norm, decreasing = TRUE))[1:1000]

set.seed(42)
boruta_res <- Boruta(x = expr_norm[, top1000], y = y, maxRuns = 200)
boruta_genes <- getSelectedAttributes(boruta_res, withTentative = FALSE)

# Union of both methods
all_selected <- unique(c(lasso_genes, boruta_genes))

# ================================================================
# STEP 7: Prepare final data with gene symbols
# ================================================================
expr_final <- expr_norm[, all_selected]
gene_map <- data.frame(ensembl = rownames(rowData(rna_data)),
                        symbol = rowData(rna_data)$Symbol)
col_symbols <- gene_map$symbol[match(colnames(expr_final), gene_map$ensembl)]
col_symbols[is.na(col_symbols)] <- colnames(expr_final)[is.na(col_symbols)]
colnames(expr_final) <- make.unique(col_symbols)

# ================================================================
# STEP 8: Discretize and run R.ROSETTA
# ================================================================
expr_disc <- as.data.frame(expr_final)
for (col in colnames(expr_disc)) {
  brks <- quantile(expr_disc[[col]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  brks[1] <- -Inf; brks[length(brks)] <- Inf
  expr_disc[[col]] <- cut(expr_disc[[col]], breaks = brks,
                           labels = c("low","medium","high"), include.lowest = TRUE)
}
expr_disc$decision <- as.factor(response_final)

set.seed(42)
result <- rosetta(expr_disc, classifier = "StandardVoter",
                   cvNum = 10, discrete = TRUE)
print(result$quality)

# ================================================================
# STEP 9: Permutation test
# ================================================================
n_perms <- 100
perm_acc <- numeric(n_perms)
for (i in 1:n_perms) {
  perm_data <- expr_disc
  perm_data$decision <- sample(perm_data$decision)
  perm_res <- tryCatch(
    rosetta(perm_data, classifier="StandardVoter", cvNum=10, discrete=TRUE),
    error = function(e) NULL)
  if (!is.null(perm_res)) perm_acc[i] <- perm_res$quality$accuracyMean
}
cat("Permutation p-value:", mean(perm_acc >= result$quality$accuracyMean), "\n")
