##### Function to find the longest run of consecutive NAs ###########

longest_na_run <- function(x) {
  # Convert to logical: TRUE if NA, FALSE otherwise
  is_na <- is.na(x)
  
  # Find lengths of all runs of TRUE/FALSE
  rle_result <- rle(is_na)
  
  # Get lengths of runs where value is TRUE (i.e., runs of NAs)
  na_lengths <- rle_result$lengths[rle_result$values == TRUE]
  
  if (length(na_lengths) == 0) {
    return(0)  # No missing values
  }
  
  # Return the maximum length
  max_na_run <- max(na_lengths)
  
  return(max_na_run)
}

