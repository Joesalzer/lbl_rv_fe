###########################
# helper functions
###########################

require(tidyverse)
require(stringr)
require(Matrix)
require(data.table)
#library(nlme)

###########################
## fit LS periodigram ##
###########################

# times : vector of time-points
# responses : vector of responses
# min_per : min period
# max_per : max period
# w : optional vector of weights
# covar_mat : optional matrix of covariates to control for
# period_num_grid : number of grid values from min to max to test the period on

ls_periodigram = function(times, response, min_per, max_per, 
                                    w = NULL, covar_mat = NULL, 
                                    grid_points = 3000,
                                    n_fourier = 1) {
  # validate n_fourier
  if (!is.numeric(n_fourier) || n_fourier < 1 || n_fourier != floor(n_fourier)) {
    stop("n_fourier must be a positive integer >= 1.")
  }
  
  # number of time-points
  n_T = length(times)
  
  # handle Weights
  if (is.null(w)) {
    w = rep(1, n_T)
  }
  # normalize weights to sum to n_T
  w = w * (n_T / sum(w)) 
  
  # setup null model 
  if (is.null(covar_mat)) {
    X0 = matrix(1, nrow = n_T, ncol = 1)
  } else {
    covar_mat = as.matrix(covar_mat)
    if (nrow(covar_mat) != n_T) stop("Covariate matrix rows must equal length of times.")
    X0 = cbind(1, covar_mat)
  }
  
  # fit null model
  null_fit = lm.wfit(x = X0, y = response, w = w)
  rss0 = sum(w * null_fit$residuals^2)
  df0  = null_fit$df.residual
  
  # log grid in period 
  periods = exp(seq(log(min_per), log(max_per), length.out = grid_points))
  freqs = 1/periods
  
  # pre-allocate storage 
  powers = numeric(grid_points)
  amps   = numeric(grid_points)
  
  # store the BEST model 
  best_rss = Inf
  best_model_coeffs = NULL
  best_period = NA
  best_lm_wfit = NULL  # <--
  
  for (i in seq_along(freqs)) {
    omega = 2 * pi * freqs[i]
    
    harmonic_cols = vector("list", 2 * n_fourier)
    for (k in seq_len(n_fourier)) {
      harmonic_cols[[2*k - 1]] = cos(k * omega * times)
      harmonic_cols[[2*k]]     = sin(k * omega * times)
    }
    harmonic_mat = do.call(cbind, harmonic_cols)
    X1 = cbind(X0, harmonic_mat)
    
    alt_fit = lm.wfit(x = X1, y = response, w = w)
    rss1 = sum(w * alt_fit$residuals^2)
    
    powers[i] = (rss0 - rss1) / rss0
    
    coefs = alt_fit$coefficients
    n_X0  = ncol(X0)
    beta_cos1 = coefs[n_X0 + 1]
    beta_sin1 = coefs[n_X0 + 2]
    amps[i] = sqrt(beta_cos1^2 + beta_sin1^2)
    
    if (rss1 < best_rss) {
      best_rss = rss1
      best_model_coeffs = coefs
      best_period = periods[i]
      best_lm_wfit = alt_fit
    }
  }
  
  return(list(
    periods   = periods,
    freqs     = freqs,
    powers    = powers,
    amps      = amps,
    n_fourier = n_fourier,
    best_fit  = list(
      period   = best_period,
      coeffs   = best_model_coeffs,
      rss      = best_rss,
      lm_wfit  = best_lm_wfit 
    )
  ))
}


###########################
## fit a LM using model matrix and responses##
###########################

# X (n x p) : model matrix
# Y (n x 1) : responses
# PRINT_TIME :logical for whether we print of not

# returns a list of model fits
sparseLM = function(X, Y, PRINT_TIME = T) {
  START_TIME = Sys.time()
  
  # crossproducts
  XtX = crossprod(X)
  XtX = Matrix::forceSymmetric(XtX, uplo = "L")
  XtY = crossprod(X, Y)
  
  # optional for sparse stability
  # eps = 1e-11
  # XtX = XtX + eps * Matrix::Diagonal(ncol(XtX))
  
  # cholesky factorization
  chol_XtX = try(Matrix::Cholesky(XtX), silent = TRUE)
  if (inherits(chol_XtX, "try-error")) {
    cat("XtX not positive definite\n")
    return(NULL)
  }
  
  # solve for beta_hat
  beta_hat = Matrix::solve(chol_XtX, XtY)
  
  # fitted values and residuals
  y_hat = X %*% beta_hat
  resid = Y - y_hat
  
  # dimensions
  p = ncol(X)
  n = nrow(X)
  
  # sum of squared residuals
  SSR = drop(crossprod(resid))
  
  # sigma2 hat and variance of beta_hat
  sigma2_hat = SSR / n
  
  # BIC calculations
  BIC = n * log(SSR/n) + (p + 1) * log(n)
  # residual standard errors
  RSE = sqrt( SSR/(n-p) )
  
  if (PRINT_TIME) {
    cat("time to fit ols model =",Sys.time() - START_TIME,"\n")
  }
  
  return( list( beta_hat = beta_hat,
                y_hat = y_hat,
                resid = resid,
                sigma2_hat = sigma2_hat,
                BIC = BIC,
                RSE = RSE))
}

###########################
## fit a weighted LM using model matrix and responses##
###########################

sparseWLM = function(X, Y, w = NULL, PRINT_TIME = TRUE) {
  START_TIME = Sys.time()
  
  # default weights
  if (is.null(w)) {
    w = rep(1, length(Y))
  }
  
  # apply sqrt weights
  sqrt_w = sqrt(w)
  
  # row-scale X and Y
  Xw = X * sqrt_w
  Yw = Y * sqrt_w

  # crossproducts
  XtWX = crossprod(Xw)
  XtWX = Matrix::forceSymmetric(XtWX, uplo = "L")
  XtWY = crossprod(Xw, Yw)

  # cholesky factorization
  chol_XtWX = try(Matrix::Cholesky(XtWX), silent = TRUE)
  if (inherits(chol_XtWX, "try-error")) {
    cat("XtWX not positive definite (even after ridge)\n")
    return(NULL)
  }

  # solve for beta_hat
  beta_hat = Matrix::solve(chol_XtWX, XtWY)

  # fitted values and residuals
  y_hat = X %*% beta_hat
  resid = Y - y_hat
  
  
  # dimensions
  p = ncol(X)
  n = nrow(X)
  
  # weight matrix
  W = Matrix::Diagonal(x = w)
  
  # effective sample size
  eff = sum(w)^2 / sum(w^2)
  
  # sigma2 hat and variance of beta_hat
  SSR = drop( t(resid) %*% W %*% resid )
  if (eff < p) {
    cat("Warning, effective sample size is smaller than number of parameters. Consider using fewer parameters. \n")
    sigma2_hat = NA
    RSE = NA
  } else {
    sigma2_hat = SSR / (eff - p)
    RSE = sqrt(sigma2_hat) 
  }
  
  # AIC/BIC calculation
  BIC = n * log(SSR / n) - sum(log(w)) + (p + 1) * log(n) + n * (1 + log(2 * pi))
  AIC = n * log(SSR / n) - sum(log(w)) + (p + 1) * 2 + n * (1 + log(2 * pi))
  
  if (PRINT_TIME) {
    print(Sys.time() - START_TIME)
  }
  
  return( list( beta_hat = beta_hat,
                y_hat = y_hat,
                resid = resid,
                sigma2_hat = sigma2_hat,
                AIC = AIC,
                BIC = BIC,
                RSE = RSE,
                eff = eff))
}

###########################
## fit IRLS model ##
###########################

# rse_tol : tolerance for the RSE as it changes
# max.iter : maximum number of iterations to run
# X (n x p) : model matrix
# Y (n x 1) : responses
# df : dataframe of observations
# lineIDname : column name of the lines

IRLS = function(X, Y, df, lineIDname, bic_tol = 1, max.iter = 10) {
  
  iter = 1
  n = nrow(X)
  p = ncol(X)
  w_irls = rep(1, n) 
  wls_fit = sparseWLM(X, Y, w_irls, PRINT_TIME = F)
  
  bic_prev = Inf
  bic_curr = wls_fit$BIC
  
  # Convert to data.table once
  dt = as.data.table(df)
  setnames(dt, lineIDname, "line_id")
  
  while( (iter < max.iter) & ( abs(bic_prev - bic_curr) > bic_tol) ) {
    cat("  Iter", iter, "| BIC:", bic_curr, "| diff BIC:", abs(bic_curr - bic_prev), "\n")
    
    iter = iter + 1

    dt[, r := wls_fit$resid[,1]]
    dt[, sar := sum(abs(r)), by = line_id]
    w_irls = 1 / pmax(dt$sar, .Machine$double.eps)
    w_irls = n * (w_irls / sum(w_irls)) # normalize to sum to n
    
    wls_fit = sparseWLM(X, Y, w_irls, PRINT_TIME = F)
    bic_prev = bic_curr
    bic_curr = wls_fit$BIC
    
  }
  
  cat("last iter: ", iter, " diff BIC: ", abs(bic_curr - bic_prev), "\n")
  return(list(fit = wls_fit, weights = w_irls))
}


###########################
# create windows to do cv over
###########################

## input ##
# x : vector of time-points
# window_size : window size of the sliding window-cv

## output ##
# list of size equal to the length of times. each index contains the test point (test_timepoint) and the set of times around that time-point to remove during the cv (test_set)
create_sliding_windows = function(x, window_size) {
  # initialize list of validations sets
  val_sets = list()
  
  # loop through all days and create windows
  for (tp in x) {
    # get the days in our dataset that are within that window
    val_set = list(test_set = x[(x >= tp - window_size) & (x <= tp + window_size)],
                   test_timepoint = as.Date(tp))
    
    val_sets = append(val_sets, list(val_set) )
  }
  
  return(val_sets)
}

###########################
## downsample a wavelength grid ##
###########################

# x : vector of wavelengths
# lower : lower limit of grid
# upper : upper limit of grid
# n : number of points to downsample

# get a filtered grid with a specific number of points between an lower and upper limit
filter_grid = function(x, lower = -Inf, upper = Inf, n = 100) {
  in_range = (x >= lower) & (x <= upper)
  idx = which(in_range)
  selected = idx[round(seq(1, length(idx), length.out = n))]
  mask = seq_along(x) %in% selected
  mask
}