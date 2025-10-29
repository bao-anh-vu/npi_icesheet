
crps_multivariate <- function(pred_samples, observed) {
  # Check dimensions
  n <- nrow(pred_samples)
  d <- ncol(pred_samples)
  
  # Ensure observed is a row vector
  if (length(observed) != d) {
    stop("Observed vector must have the same dimensions as the prediction samples.")
  }
  
  # Calculate the first term: Mean L2 norm between each sample vector and observed vector
  term1 <- mean(apply(pred_samples, 1, function(x) sqrt(sum((x - observed)^2))))
  
  # Calculate the second term: Mean L2 norm between all pairs of sample vectors
  term2 <- mean(apply(pred_samples, 1, function(x1) {
    mean(apply(pred_samples, 1, function(x2) sqrt(sum((x1 - x2)^2))))
  }))
  
  # Energy Score (Multivariate CRPS) is the difference between these two terms
  energy_score <- term1 - (term2 / 2)
  return(energy_score)
}

## Compute RMSE
rmse <- function(obs, sim) {
  sqrt(mean((obs - sim)^2, na.rm = T))
}
