# ===================================================================
# Function: compute_pooled_mase()
# Purpose : For each imputation and each shift q (-k_t to +k_t),
#           computes the MASE between the imputed target series
#           and the shifted S0 vector, then pools the results
#           across the m imputations.
#
# Inputs:
#   Pr_imputed : matrix/data.frame with m columns (the m imputations of your target)
#   shift      : matrix with (2*k_t + 1) columns (the shifted S0 vectors)
#   m          : number of imputations (default: ncol(Pr_imputed))
#   k_t        : maximum shift (default: derived from ncol(shift))
#   step       : step size for MASE (e.g. 1 or 24)
#
# Returns: a list containing
#   $MASE        : matrix (m rows × (2*k_t+1) columns)
#   $PooledMASE  : named vector of length (2*k_t+1) with the averaged MASE per shift
# ===================================================================

compute_pooled_mase <- function(Pr_imputed, 
                                shift, 
                                step, 
                                m = NULL, 
                                k_t = NULL) {
  
  # Derive m and k_t if not supplied
  if (is.null(m))   m   <- ncol(Pr_imputed)
  if (is.null(k_t)) k_t <- (ncol(shift) - 1) / 2
  
  n_shifts <- 2 * k_t + 1
  
  # Initialise MASE matrix
  MASE <- matrix(NA, nrow = m, ncol = n_shifts)
  colnames(MASE) <- paste0("Shift_", (-k_t):k_t)
  
  # Compute MASE for every imputation and every shift
  for (imp in 1:m) {
    for (s in 1:n_shifts) {
      df_temp <- na.omit(cbind(Pr_imputed, shift[, s]))
      
      # actual = shifted S0 (last column), predicted = current imputation
      MASE[imp, s] <- mase(
        actual    = df_temp[, ncol(df_temp)],
        predicted = df_temp[, imp],
        step      = step
      )
    }
  }
  
  # Pool across imputations
  PooledMASE <- colMeans(MASE)
  
  return(list(
    MASE        = MASE,
    PooledMASE  = PooledMASE
  ))
}
