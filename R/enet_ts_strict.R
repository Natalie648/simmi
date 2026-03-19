############### Function to define custom elastic net method for MICE using time-series CV #######

# Using time-series block CV for BOTH hyper-parameters
# Clean joint hyperparameter optimization.
# True joint tuning involves:
# 1. Define block folds
# For each fold:
#  For each α: Fit glmnet over λ path
#  Evaluate prediction error on that fold
# Aggregate error across folds
# Select best (α, λ) pair globally
# Requires writing manual CV logic rather than relying on cv.glmnet().
# Written for stochastic imputation without bootstrapping

mice.impute.enet_ts_strict <- function(y, ry, x, ...) {
  
  y_obs <- y[ry]
  x_obs <- as.matrix(x[ry, , drop = FALSE])
  x_mis <- as.matrix(x[!ry, , drop = FALSE])
  
  n_obs <- length(y_obs)
  
  foldid <- make_time_folds(n_obs)
  n_splits=nfolds
  
  alpha_grid <- c(0.1, 0.5, 0.9, 1.0)
  
  best_error <- Inf
  best_alpha <- NULL
  best_lambda <- NULL
  
  for (a in alpha_grid) {
    
    # Fit on full to get common lambda path
    fit_full <- glmnet(x_obs, y_obs, alpha = a)
    lambda_path <- fit_full$lambda
    n_lambda <- length(lambda_path)
    
    cv_errors <- matrix(NA, nrow = n_splits, ncol = n_lambda)
    
    for (k in 1:n_splits) {
      test_idx <- which(foldid == k)
      train_idx <- which(foldid != k)
      
      if (length(train_idx) < ncol(x_obs) + 1 || length(test_idx) == 0) next  # Skip if train too small or test empty
      
      # Fit on train with same lambda path
      fit_k <- glmnet(x_obs[train_idx, ], y_obs[train_idx], alpha = a, lambda = lambda_path)
      
      preds <- predict(fit_k, newx = x_obs[test_idx, ])
      
      mse <- colMeans((y_obs[test_idx] - preds)^2)
      
      cv_errors[k, ] <- mse
    }
    
    mean_mse <- colMeans(cv_errors, na.rm = TRUE)
    
    if (!all(is.na(mean_mse))) {
      min_mse <- min(mean_mse, na.rm = TRUE)
      
      if (min_mse < best_error) {
        best_error <- min_mse
        best_alpha <- a
        best_lambda <- lambda_path[which.min(mean_mse)]
      }
    }
  }
  
  if (is.null(best_alpha)) {
    # Fallback to lasso with minimal regularization
    best_alpha <- 1.0
    fallback_fit <- glmnet(x_obs, y_obs, alpha = best_alpha)
    best_lambda <- fallback_fit$lambda[length(fallback_fit$lambda)]  # Smallest lambda
  }
  
  final_fit <- glmnet(x_obs, y_obs, alpha = best_alpha, lambda = best_lambda)
  
  # Deterministic predictions
  y_pred_det <- predict(final_fit, newx = x_mis)[, 1]
  
  # For stochastic imputation: add noise from residual distribution
  preds_obs <- predict(final_fit, newx = x_obs)[, 1]
  res <- y_obs - preds_obs
  df_eff <- final_fit$df + 1  # Effective df including intercept
  denom <- n_obs - df_eff
  sigma_hat <- if (denom > 0) {
    sqrt(sum(res^2) / denom)
  } else {
    sd(res)  # Fallback
  }
  
  epsilon <- rnorm(nrow(x_mis), mean = 0, sd = sigma_hat)
  y_pred <- y_pred_det + epsilon
  
  return(as.numeric(y_pred))
}
