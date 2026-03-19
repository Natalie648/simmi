# ===================================================================
# Function: generate_pr_imputed()
# Purpose : 
#   - Adds lagged/leading versions (p(r) vectors) of all numeric predictors
#   - Subsets to the spatial window [l_start : l_end]
#   - Keeps only columns with ≥ 40% observed values
#   - Blanks out the irregular pattern period in the target column (column 1)
#   - Runs MICE using custom time-series method ("enet_ts_strict" or "lasso_ts")
#   - Returns the m imputed versions of the target series
#
# Inputs:
#   my_spatial_ts : cleaned and prepared data frame (target in column 1, factors already created)
#   l_start, l_end: row indices of the window
#   l             : exact window length (from compute_spatial_window)
#   k             : length of the irregular pattern
#   p             : lag/lead range (e.g. 2 → shifts -2:2)
#   m             : number of imputations
#   method        : "enet_ts_strict" or "lasso_ts"
#   seed          : random seed for reproducibility (default 12345)
#
# Returns: a data frame/matrix with m columns = the m imputed versions of the target series
# ===================================================================

generate_pr_imputed <- function(my_spatial_ts, 
                                l_start, l_end, l, k, p, m, 
                                method, 
                                seed = 12345) {
  
  # Work on a local copy (never modifies the original)
  df <- my_spatial_ts
  
  # -----------------------------------------------------------------
  # 1. Build p(r) vectors: keep original target + shifted numeric predictors + factors
  # -----------------------------------------------------------------
  num_numeric <- sum(sapply(df, is.numeric))
  
  shifted_predictors <- data.table::shift(
    df[, c(2:num_numeric), drop = FALSE],
    n = -p:p,
    type = "shift",
    give.names = TRUE
  )
  
  non_numeric <- df[, !sapply(df, is.numeric), drop = FALSE]
  
  df <- cbind(
    df[, 1, drop = FALSE],          # target column (p0)
    shifted_predictors,
    non_numeric
  )
  
  # -----------------------------------------------------------------
  # 2. Subset to the spatial window
  # -----------------------------------------------------------------
  df <- df[l_start:l_end, , drop = FALSE]
  
  # -----------------------------------------------------------------
  # 3. Keep only columns with at least 40% observed values
  # -----------------------------------------------------------------
  df <- df[, colMeans(!is.na(df)) >= 0.4, drop = FALSE]
  
  # -----------------------------------------------------------------
  # 4. Blank out the irregular pattern in the target (column 1)
  # -----------------------------------------------------------------
  pad      <- (l - k) / 2
  na_start <- pad + 1
  na_end   <- pad + k
  df[na_start:na_end, 1] <- NA
  
  # -----------------------------------------------------------------
  # 5. Define imputation methods and apply custom TS method to all numeric columns
  # -----------------------------------------------------------------
  meth <- make.method(df)
  num_numeric_final <- sum(sapply(df, is.numeric))
  meth[1:num_numeric_final] <- method
  
  # -----------------------------------------------------------------
  # 6. Run MICE
  # -----------------------------------------------------------------
  imp <- mice(
    data   = df,
    m      = m,
    maxit  = 10,
    method = meth,
    seed   = seed,
    print  = FALSE
  )
  
  # -----------------------------------------------------------------
  # 7. Extract the m imputations of the target series
  # -----------------------------------------------------------------
  data_imp   <- complete(imp, "repeated")
  Pr_imputed <- data_imp[, 1:m, drop = FALSE]
  
  return(Pr_imputed)
}