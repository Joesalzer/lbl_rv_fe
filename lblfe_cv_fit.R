# ── pointer to source of helper ───────────────────────────────────────────────
suppressPackageStartupMessages({
  ## local ##
  source("rv_helper.R")
  
  ## CHTC ##
  #source(paste0(getwd(),"/rv_helper.R"))
})

# ── command-line arguments ─────────────────────────────────────────────────────
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript lblfe_cv_fit.R <inputCSVPath> <outputDir>\n",
       "Example: Rscript lblfe_cv_fit.R /path/to/input.csv /path/to/output/dir")
}
inputCSVPath = args[1]
outputDir    = args[2]

# ── directories & file paths ───────────────────────────────────────────────────

# create output directory if it does not exist
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

fileName     = basename(inputCSVPath)
results_file = sub("\\.csv$", "", file.path(outputDir, fileName))

# ── read in data ───────────────────────────────────────────────────────────────
lbl_df = fread(inputCSVPath) %>%
  mutate(lineID = factor(lineID),
         timeID = factor(timeID),
         date   = as.Date(date))

# ── validate & extract column names ───────────────────────────────────────────
expected_cols   = c("lineID", "timeID", "date", "split", "rv_dirty")
lbl_df_colnames = colnames(lbl_df)

if (!identical(lbl_df_colnames[1:5], expected_cols)) {
  stop(
    "First five columns are not correctly specified.\n",
    "Expected: ", paste(expected_cols,        collapse = ", "), "\n",
    "Got:      ", paste(lbl_df_colnames[1:5], collapse = ", ")
  )
}

# save remaining column names after the first five (empty character vector if none)
covar_colnames = lbl_df_colnames[seq_along(lbl_df_colnames) > 5]

# ── split into train and test_timepoint ───────────────────────────────────────

# training rows: fit the model on these
train_df = lbl_df %>%
  filter(split == "train") %>%
  mutate(lineID = factor(lineID),
         timeID = factor(timeID),
         date   = as.Date(date))

# test timepoint rows: evaluate on these (one observation per line)
test_df  = lbl_df %>%
  filter(split == "test_timepoint") %>%
  mutate(lineID = factor(lineID),
         timeID = factor(timeID),
         date   = as.Date(date))

cat("Training rows    :", nrow(train_df), "\n")
cat("Test TP     :", as.character(test_df$date[1]),  "\n")
cat("Test TP rows     :", nrow(test_df),  "\n")
cat("Covariates       :", if (length(covar_colnames) == 0) "none" else paste(covar_colnames, collapse = ", "), "\n\n")

# ── build training design matrix ──────────────────────────────────────────────
if (length(covar_colnames) == 0) {
  formula_terms = c("timeID", "lineID")
} else {
  formula_terms = c("timeID", "lineID", paste("lineID", covar_colnames, sep = ":"))
}

model_formula = reformulate(formula_terms)

train_designMat = sparse.model.matrix(
  model_formula,
  train_df,
  contrasts.arg = list(timeID = "contr.sum", lineID = "contr.sum")
)

cat("Training design matrix dimensions:", dim(train_designMat), "\n\n")

# ── auxiliary objects ──────────────────────────────────────────────────────────
train_obsIndices = train_df %>%
  dplyr::select(lineID, timeID)

train_responses = train_df$rv_dirty
lineIDs         = levels(lbl_df$lineID)
n_L             = length(lineIDs)

# ── helper: predict at test_timepoint ─────────────────────────────────────────
# since test_timepoint has no time fixed effect, we build a design matrix
# using only lineID and covariate terms (no timeID) and predict from those
predict_test_tp = function(beta_hat, test_df, covar_colnames) {
  
  # build test design matrix with only lineID and covariate terms (no timeID)
  if (length(covar_colnames) == 0) {
    test_formula_terms = c("lineID")
  } else {
    test_formula_terms = c("lineID", paste("lineID", covar_colnames, sep = ":"))
  }
  
  test_formula = reformulate(test_formula_terms)
  
  # use only lineID levels seen in training to avoid new factor level issues
  test_designMat = sparse.model.matrix(
    test_formula,
    test_df,
    contrasts.arg = list(lineID = "contr.sum")
  )
  
  # extract beta coefficients corresponding to columns in test design matrix
  # time fixed effect columns in beta_hat are excluded since they are not in test_designMat
  shared_cols  = intersect(colnames(test_designMat), rownames(beta_hat))
  beta_test    = beta_hat[shared_cols, , drop = FALSE]
  test_designMat_shared = test_designMat[, shared_cols, drop = FALSE]
  
  # predicted rv at test timepoint (line FE + covariate terms only)
  rv_pred = as.numeric(test_designMat_shared %*% beta_test)
  
  return(rv_pred)
}

# ── OLS ────────────────────────────────────────────────────────────────────────
cat("Fitting OLS model...\n")
ols_fit  = sparseWLM(train_designMat, train_responses, w = NULL, PRINT_TIME = TRUE)
#ols_fit = sparseLM(train_designMat, train_responses, PRINT_TIME = TRUE)

cat("## OLS RESULTS ##\n",
    str_c(" BIC: ",        formatC(ols_fit$BIC)),    "\n",
    str_c(" RSE: ",        round(ols_fit$RSE, 4)),   "\n",
    str_c(" Parameters: ", ncol(train_designMat)),   "\n",
    str_c(" Sample size: ", nrow(train_designMat)),  "\n\n")

# predict at test timepoint using line FE and covariates only
ols_rv_pred = predict_test_tp(ols_fit$beta_hat, test_df, covar_colnames)

# build results data.table with lineID, timeID, date, rv_dirty, rv_pred, resid
ols_results_dt = test_df %>%
  dplyr::select(lineID, timeID, date, rv_dirty) %>%
  mutate(rv_pred = ols_rv_pred,
         resid   = rv_dirty - rv_pred)

saveRDS(list(fit        = ols_fit,
             obsIndices = train_obsIndices),
        str_c(results_file, "__fit=OLS.rds"))
fwrite(ols_results_dt, str_c(results_file, "__fit=OLS__predictions.csv"))

# ── IRLS ───────────────────────────────────────────────────────────────────────
cat("Fitting IRLS model...\n")
irls_results = IRLS(bic_tol = 1, max.iter = 45,
                    X = train_designMat, Y = train_responses,
                    df = train_obsIndices, lineIDname = "lineID")
irls_fit = irls_results$fit
w_irls   = irls_results$weights

cat("## IRLS RESULTS ##\n",
    str_c(" BIC: ",                   formatC(irls_fit$BIC)),        "\n",
    str_c(" RSE: ",                   round(irls_fit$RSE, 4)),       "\n",
    str_c(" Parameters: ",            ncol(train_designMat)),          "\n",
    str_c(" Effective sample size: ", round(1 / sum(w_irls^2), 1)),    "\n\n")

# predict at test timepoint using line FE and covariates only
irls_rv_pred = predict_test_tp(irls_fit$beta_hat, test_df, covar_colnames)

# build results data.table with lineID, timeID, date, rv_dirty, rv_pred, resid
irls_results_dt = test_df %>%
  dplyr::select(lineID, timeID, date, rv_dirty) %>%
  mutate(rv_pred = irls_rv_pred,
         resid   = rv_dirty - rv_pred)

saveRDS(list(fit        = irls_fit,
             weights    = w_irls,
             obsIndices = train_obsIndices),
        str_c(results_file, "__fit=IRLS.rds"))
fwrite(irls_results_dt, str_c(results_file, "__fit=IRLS__predictions.csv"))