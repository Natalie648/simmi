# ===================================================================
# Function: prep_spatial_ts()
# Purpose : Takes the original data frame (my_spatial_ts_orig) and
#           returns a prepared version (my_spatial_ts) ready for MICE.
#           - Removes all Date / datetime columns
#           - Removes any logical columns 
#           - Does not overwrite the original data frame
# ===================================================================

prep_spatial_ts <- function(my_spatial_ts_orig) {
  
  # Step 1: Drop Date and datetime columns
  my_spatial_ts <- my_spatial_ts_orig %>%
    select(!where(is.Date) & !where(lubridate::is.instant))
  
  # Step 2: Drop logical columns
  my_spatial_ts <- my_spatial_ts[, !(sapply(my_spatial_ts, is.logical))]
  
  return(my_spatial_ts)
}