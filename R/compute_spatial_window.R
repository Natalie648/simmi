# ===================================================================
# Function: compute_spatial_window()
# Purpose : Calculates the optimal window length `l` and the bounds
#           `l_start` / `l_end` so that the irregular pattern
#           (from k_start to k_end) is nicely centered.
#           Automatically clips to the data boundaries and shrinks
#           k_t when the window is too small for the full shifts.
#
# Inputs  :
#   my_spatial_ts : your cleaned and prepared data frame
#   k_t           : desired shift (e.g. 5*theta)
#   k_start       : row number where the irregular pattern STARTS
#   k_end         : row number where the irregular pattern ENDS
#
# Returns : a list with the final values of l, l_start, l_end, k_t
# ===================================================================

compute_spatial_window <- function(my_spatial_ts, k_t, k_start, k_end) {
  
  k <- k_end-k_start+1
  n_rows <- nrow(my_spatial_ts)
  
  # 1. Choose window length l 
  l <- ifelse(((ncol(my_spatial_ts) - 1) * 20) > 2 * k_t,
              k + ((ncol(my_spatial_ts) - 1) * 20),
              k + 2 * k_t)
  
  # 2. Center the window around the irregular pattern
  l_start <- k_start - ((l - k) / 2)
  l_end   <- k_end   + ((l - k) / 2)
  
  # 3. Clip start to row 1
  l_start <- ifelse(l_start < 1, 1, l_start)
  
  # 4. Tentative end based on new start
  l_end <- l + l_start - 1
  
  # 5. Clip end to the last row of the data
  l_end <- ifelse(l_end > n_rows, n_rows, l_end)
  
  # 6. Adjust start again so the window still has length l
  l_start <- max(l_end - l + 1, 1)
  
  # 7. Final actual window length
  l <- l_end - l_start + 1
  
  # 8. Reduce k_t if the window is too short for full shifts
  k_t <- min(k_t, k_start - l_start, l_end - k_end)
  
  # Return everything in a list
  return(list(
    k       = k,
    l       = l,
    l_start = l_start,
    l_end   = l_end,
    k_t     = k_t
  ))
}
