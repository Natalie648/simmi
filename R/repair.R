# ===================================================================
# Function: repair_spatial_ts()
# Purpose : Creates a copy of the original data frame and replaces
#           ONLY the irregular pattern period (k_start:k_end) in
#           column d with the optimally shifted version of itself.
#           Uses data.table::shift directly on the original series.
#
# Inputs:
#   my_spatial_ts_orig : original full data frame (unchanged)
#   k_start            : starting row of the irregular pattern
#   k_end              : ending row of the irregular pattern
#   d                  : column index or name of the target series
#   shifts             : vector of possible shifts (-k_t:k_t)
#   opt_shift_idx      : index of the best shift (e.g. from which.min(PooledMASE))
#
# Returns: my_spatial_ts_repaired = full original data frame with
#          the irregular period now filled with the optimal shift
# ===================================================================

repair_spatial_ts <- function(my_spatial_ts_orig,
                              k_start,
                              k_end,
                              d,
                              shifts,
                              opt_shift_idx) {
  
  # Safety checks
  if (k_start < 1 || k_end > nrow(my_spatial_ts_orig)) {
    stop("k_start:k_end is outside the rows of my_spatial_ts_orig")
  }
  if (opt_shift_idx < 1 || opt_shift_idx > length(shifts)) {
    stop("opt_shift_idx is invalid for the shifts vector")
  }
  
  # 1. Extract the optimal shift amount
  shift_amount <- shifts[opt_shift_idx]
  
  # 2. Create the replacement vector (exactly length k)
  replacement_vec <- as.numeric(
    data.table::shift(
      my_spatial_ts_orig[k_start:k_end, d, drop = TRUE],
      n          = shift_amount,
      type       = "shift",
      give.names = FALSE
    )
  )
  
  # 3. Create repaired copy and insert the vector
  my_spatial_ts_repaired <- my_spatial_ts_orig
  my_spatial_ts_repaired[k_start:k_end, d] <- replacement_vec
  
  return(my_spatial_ts_repaired)
}