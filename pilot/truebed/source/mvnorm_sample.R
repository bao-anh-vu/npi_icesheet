#Test: sampling from a multivariate normal
mvnorm_sample <- function(n, m, V) {
  samples <- matrix(NA, nrow = nrow(V), ncol = n)

  if (all(abs(V) < .Machine$double.eps)) { # Check if V is all zero
    samples <- matrix(rep(m, n), nrow = nrow(V), ncol = n) # just repeat the mean n times
  } else {
    L <- chol(V)
    for (i in 1:n) {
      samples[, i] <- m + L %*% rnorm(nrow(V))
    }
  }
  samples
}
