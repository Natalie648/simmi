# ===================================================================
# Function: compute_mase_significance()
# Purpose : Performs one-sided significance testing on the Pooled MASE
#           values using Newey-West HAC standard errors. Applies Rubin's 
#           rules to combine within- and between-imputation variance, then 
#           tests H0: MASE = null vs HA: MASE < null.
#
# Inputs:
#   Pr_imputed   : matrix with m columns (the m imputations of the target)
#   shift        : matrix with (2*k_t + 1) columns (shifted S0 vectors)
#   MASE         : matrix (m × (2*k_t+1)) from compute_pooled_mase()
#   PooledMASE   : named vector of pooled MASE values
#   m            : number of imputations
#   k_t          : maximum shift
#   step         : step size used in MASE
#   null         : null value for MASE (default = 1)

#
# Returns: a list containing
#   $w_vec     : within-imputation variance (for MAE)
#   $b_vec     : between-imputation variance (for MAE)
#   $T_MAE     : total variance for MAE_bar (Rubin's rules)
#   $T_MASE    : total variance for MASE_bar
#   $df_vec    : Rubin’s degrees of freedom
#   $stat_vec  : t-statistics
#   $p_vec     : one-sided p-values (lower tail)
# ===================================================================

compute_mase_significance <- function(Pr_imputed,
                                      shift,
                                      MASE,
                                      PooledMASE,
                                      m,
                                      k_t,
                                      step,
                                      null) {
  
  n_shifts <- 2 * k_t + 1
  
  w_vec     <- numeric(n_shifts)   # within-var for MAE
  b_vec     <- numeric(n_shifts)   # between-var for MAE
  T_MAE_vec <- numeric(n_shifts)   # total var for MAE_bar
  T_MASE_vec<- numeric(n_shifts)
  df_vec    <- numeric(n_shifts)
  stat_vec  <- numeric(n_shifts)
  p_vec     <- numeric(n_shifts)
  
  for (i in 1:n_shifts) {
    overlap_rows <- which(!is.na(shift[, i]))
    n <- length(overlap_rows)
    
    if (n < 50) {                     
      p_vec[i] <- NA
      next
    }
    
    actual <- shift[overlap_rows, i]
    naive  <- mean(abs(diff(actual, lag = step)))   # denominator of MASE
    
    if (naive == 0 || is.na(naive)) {
      p_vec[i] <- NA
      next
    }
    
    # Back-compute MAE_k from MASE_k * naive
    MAE <- MASE[, i] * naive
    
    # Within-imputation variance using Newey-West HAC on each imputation
    var_a <- numeric(m)
    for (a in 1:m) {
      pred       <- Pr_imputed[overlap_rows, a]
      errors     <- actual - pred
      abs_errors <- abs(errors)
      
      # HAC variance of the mean absolute error
      lm_mod <- lm(abs_errors ~ 1)
      hac_var <- sandwich::NeweyWest(
        lm_mod,
        lag      = min(floor(bwNeweyWest(lm_mod)),floor(4*(n/100)^(2/9))),      # ←←← Adjustable via function argument
        prewhite = FALSE,
        adjust   = TRUE
      )[1, 1]
      
      var_a[a] <- hac_var
    }
    
    w_vec[i]     <- mean(var_a, na.rm = TRUE)
    b_vec[i]     <- var(MAE, na.rm = TRUE)
    T_MAE        <- w_vec[i] + (1 + 1/m) * b_vec[i]
    T_MAE_vec[i] <- T_MAE                    
    T_MASE_vec[i] <- T_MAE / naive^2         
    
    if (T_MAE <= 0 || is.na(T_MAE)) {
      p_vec[i] <- NA
      next
    }
    
    T_MASE       <- T_MAE / naive^2          # Propagate to MASE scale
    
    # Rubin’s degrees of freedom
    r_i <- (1 + 1/m) * b_vec[i] / w_vec[i]
    df_vec[i] <- if (r_i <= 0 || is.infinite(r_i) || is.na(r_i)) {
      Inf
    } else {
      (m - 1) * (1 + 1/r_i)^2
    }
    
    stat_vec[i] <- (PooledMASE[i] - null) / sqrt(T_MASE)
    p_vec[i]    <- pt(stat_vec[i], df_vec[i])   # one-sided lower tail
  }
  
  return(list(
    w_vec  = w_vec,
    b_vec  = b_vec,
    T_MAE  = T_MAE_vec,
    T_MASE = T_MASE_vec,
    df_vec = df_vec,
    stat_vec = stat_vec,
    p_vec  = p_vec
  ))
}


