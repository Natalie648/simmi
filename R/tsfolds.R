############## Function to create contiguous time-series block fold IDs #######

# Ensures remainder is handled appropriately to fully maintaing time ordering
# within blocks

make_time_folds <- function(n, n_splits = nfolds) {
  fold_sizes <- rep(floor(n / n_splits), n_splits)
  remainder <- n %% n_splits
  if (remainder > 0) fold_sizes[1:remainder] <- fold_sizes[1:remainder] + 1
  
  foldid <- rep(1:n_splits, times = fold_sizes)
  return(foldid)
}
