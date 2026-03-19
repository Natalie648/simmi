##### Function to define custom LASSO method for MICE using time-series CV ####

# For stochastic imputation without bootstrapping

mice.impute.lasso_ts <- function(y, ry, x, ...) {
  
  # Observed rows
  y_obs <- y[ry]
  x_obs <- as.matrix(x[ry, , drop = FALSE])
  x_mis <- as.matrix(x[!ry, , drop = FALSE])
  
  n_obs <- length(y_obs)
  
  # Create 10 time-based folds
  foldid <- make_time_folds(n_obs, n_splits = nfolds)
  
  # Cross-validated LASSO
  cvfit <- cv.glmnet(
    x = x_obs,
    y = y_obs,
    alpha = 1, # LASSO
    foldid = foldid,
    standardize = TRUE
  )
  
  # Selected lambda
  best_lambda <- cvfit$lambda.min
  
  # Fit final model on all observed data
  final_fit <- glmnet(
    x = x_obs,
    y = y_obs,
    alpha = 1,
    lambda = best_lambda,
    standardize = TRUE
  )
  
  # Deterministic predictions
  y_pred_det <- predict(final_fit, newx = x_mis)[, 1]
  
  # For stochastic imputation: add noise from residual distribution
  preds_obs <- predict(final_fit, newx = x_obs)[, 1]
  res <- y_obs - preds_obs
  df_eff <- final_fit$df + 1 # Effective df including intercept
  denom <- n_obs - df_eff
  sigma_hat <- if (denom > 0) {
    sqrt(sum(res^2) / denom)
  } else {
    sd(res) # Fallback
  }
  
  epsilon <- rnorm(nrow(x_mis), mean = 0, sd = sigma_hat)
  y_pred <- y_pred_det + epsilon
  
  return(as.numeric(y_pred))
}
