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
  stop("Usage: Rscript lblfe_fit.R <inputCSVPath> <outputDir>\n",
       "Example: Rscript lblfe_fit.R /path/to/input.csv /path/to/output/dir")
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
expected_cols  = c("lineID", "timeID", "date", "rv_dirty")
lbl_df_colnames = colnames(lbl_df)

if (!identical(lbl_df_colnames[1:4], expected_cols)) {
  stop(
    "First four elements are not correctly specified.\n",
    "Expected: ", paste(expected_cols,        collapse = ", "), "\n",
    "Got:      ", paste(lbl_df_colnames[1:4], collapse = ", ")
  )
}
# save remaining column names after the first four (empty character vector if none)
covar_colnames = lbl_df_colnames[seq_along(lbl_df_colnames) > 4]

# ── build design matrix ────────────────────────────────────────────────────────
if (length(covar_colnames) == 0) {
  formula_terms = c("timeID", "lineID")
} else {
  formula_terms = c("timeID", "lineID", paste("lineID", covar_colnames, sep = ":"))
}

model_formula = reformulate(formula_terms)

designMat = sparse.model.matrix(
  model_formula,
  lbl_df,
  contrasts.arg = list(timeID = "contr.sum", lineID = "contr.sum")
)

cat("File name:", fileName, "\n")
cat("Design matrix dimensions:", dim(designMat), "\n")

# ── auxiliary objects ──────────────────────────────────────────────────────────
obsIndices = lbl_df %>%
  mutate(timeID = factor(timeID),
         lineID = factor(lineID)) %>%
  dplyr::select(lineID, timeID)

responses   = lbl_df$rv_dirty
time_points = unique(as.Date(lbl_df$date))
n_T         = length(time_points)
lineIDs     = unique(lbl_df$lineID)
n_L         = length(lineIDs)

# ── linear operator matrix ─────────────────────────────────────────────────────
linear_op_mat = Matrix(0, nrow = n_T, ncol = ncol(designMat), sparse = TRUE)
linear_op_mat[, 2:n_T] = contr.sum(n_T)

# ── OLS ────────────────────────────────────────────────────────────────────────
cat("Start fitting OLS model...\n")
ols_fit  = sparseWLM(designMat, responses, w = NULL, PRINT_TIME = TRUE)
#ols_fit  = sparseLM(designMat, responses, PRINT_TIME = TRUE)
cleanRV_ols  = (linear_op_mat %*% ols_fit$beta_hat[, 1])[, 1]
rms_ols      = sqrt(mean(cleanRV_ols^2))

cat("## OLS RESULTS ##\n",
    str_c(" BIC: ",        formatC(ols_fit$BIC)),          "\n",
    str_c(" RSE: ",        round(ols_fit$RSE, 4)),         "\n",
    str_c(" Parameters: ", ncol(designMat)),                   "\n",
    str_c(" Sample size: ", nrow(designMat)),                  "\n",
    str_c(" RMS: ",  round(rms_ols, 4)),                 "\n\n")

saveRDS(
  list(fit = ols_fit,
       cleanRV     = cleanRV_ols,
       time_points = time_points,
       obsIndices  = obsIndices),
  str_c(results_file,"__fit=OLS.rds")
)

# ── IRLS ───────────────────────────────────────────────────────────────────────
cat("Start fitting IRLS model...\n")
irls_results     = IRLS(bic_tol = 1, max.iter = 45, X = designMat, Y = responses, df = obsIndices, lineIDname = "lineID")
irls_fit = irls_results$fit
w_irls       = irls_results$weights
cleanRV_irls = (linear_op_mat %*% irls_fit$beta_hat[, 1])[, 1]
rms_irls     = sqrt(mean(cleanRV_irls^2))

cat("## IRLS RESULTS ##\n",
    str_c(" BIC: ", formatC(irls_fit$BIC)), "\n",
    str_c(" RSE: ", round(irls_fit$RSE, 4)), "\n",
    str_c(" Parameters: ",            ncol(designMat)), "\n",
    str_c(" Effective sample size: ", round(irls_fit$eff, 1)), "\n",
    str_c(" RMS: ",            round(rms_irls, 4)), "\n\n")

saveRDS(
  list(fit = irls_fit,
       cleanRV     = cleanRV_irls,
       time_points = time_points,
       weights     = w_irls,
       obsIndices  = obsIndices),
  str_c(results_file,"__fit=IRLS.rds")
)