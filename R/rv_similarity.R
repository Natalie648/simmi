# ===================================================================
# Function: compute_rv_similarity()
# Purpose : Calculates the RV coefficient between two multivariate data sets
#           before and after repair of the irregular pattern.
#           First multivariate data set is the one affected by irregular pattern.
#           Second multivariate data set relates to a chosen reference location.
#           RV coefficient is calculated for the period spanning the irregular 
#           pattern (from k_start to k_end).
#           Function automatically omits any rows containing NA values.
#           Function calculates RV coefficients for both univariate and multivariate data 
#
# Inputs  :
#   target_orig     : original multivariate data set containing period of irregular pattern
#   target_repaired : multivariate data set 1 after repair
#   ref_orig        : multivariate data set for reference location
#   j               : number of variables common to both data sets
#   k_start         : row number where the irregular pattern STARTS
#   k_end           : row number where the irregular pattern ENDS
#
# Returns : a table with computed RV coefficients
# ===================================================================

compute_rv_similarity <- function(target_orig, target_repaired, ref_orig, j, k_start, k_end) {
  
  # Helper: compute multivariate RV
  multivariate_rv <- function(target, ref, j, k_start, k_end) {
    comp <- cbind(target, ref)
    comp <- comp[k_start:k_end, ]
    comp <- na.omit(comp)
    similarity <- coeffRV(comp[, 1:j], comp[, (j + 1):(2*j)])
    similarity$rv
  }
  
  # Helper: compute univariate RV for each variable
  univariate_rv <- function(target, ref, j, k_start, k_end) {
    comp <- cbind(target, ref)
    comp <- comp[k_start:k_end, ]
    comp <- na.omit(comp)
    rv_values <- numeric(j)
    for(i in 1:j){
      rv_values[i] <- cor(comp[, i], comp[, j + i])^2
    }
    names(rv_values) <- paste("Variable", 1:j)
    rv_values
  }
  
  # Compute RVs
  similarity_pre <- multivariate_rv(target_orig, ref_orig, j, k_start, k_end)
  similarity_post <- multivariate_rv(target_repaired, ref_orig, j, k_start, k_end)
  
  similarity_univariate_pre <- univariate_rv(target_orig, ref_orig, j, k_start, k_end)
  similarity_univariate_post <- univariate_rv(target_repaired, ref_orig, j, k_start, k_end)
  
  # Combine univariate RVs
  similarity_univariate <- cbind(similarity_univariate_pre, similarity_univariate_post)
  colnames(similarity_univariate) <- c("RV_pre","RV_post")
  
  # Combine total RVs
  similarity_total <- cbind(similarity_pre, similarity_post)
  rownames(similarity_total) <- "Total"
  
  # Final combined table
  similarity <- rbind(similarity_univariate, similarity_total)
  
  return(similarity)
}