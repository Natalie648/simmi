# ===================================================================
# Function: generate_s0_shifts()
# Purpose : Creates the S0 vector (Series_irregular) with the irregular
#           pattern centered and padded with NAs, then generates the
#           full shift matrix s(q) for all lags from -k_t to +k_t.
#
# Inputs  :
#   my_spatial_ts : cleaned data frame (irregular series must be in column 1)
#   l_start       : starting row of the window
#   l_end         : ending row of the window
#   l             : total window length
#   k             : length of the irregular pattern
#   k_t           : maximum number of shifts in each direction
#   theta         : periodicity of the data
#
# Returns : a list containing:
#   $Series_irregular : numeric vector of length l (S0 with NAs)
#   $shift            : numeric matrix with (2*k_t + 1) columns
# ===================================================================

generate_s0_shifts <- function(my_spatial_ts, l_start, l_end, l, k, k_t, theta) {
  
  # Extract the window around the irregular pattern (column 1)
  Series_irregular <- my_spatial_ts[c(l_start:l_end), 1, drop = TRUE]
  Series_irregular <- as.numeric(Series_irregular)
  
  # Center the irregular pattern by adding NA padding on both sides
  pad <- (l - k) / 2
  Series_irregular[1:pad] <- NA
  Series_irregular[(pad + k + 1):l] <- NA
  
  # Interpolate any missing values within the irregular sequence 
  Series_irregular[(pad+1):(pad+k)] <- as.numeric(forecast::na.interp(ts(Series_irregular[(pad+1):(pad+k)], frequency=theta)))
  
  
  # Create the full shift matrix s(q): each column is a shifted version
  shift <- as.matrix(
    data.frame(
      data.table::shift(
        Series_irregular,
        n = -k_t:k_t,
        type = "shift",
        give.names = TRUE
      )
    )
  )
  
  return(list(
    Series_irregular = Series_irregular,
    shift            = shift
  ))
}