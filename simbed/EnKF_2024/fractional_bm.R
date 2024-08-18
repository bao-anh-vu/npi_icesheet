## Fractional Brownian motion ##

sim_frac_bm <- function(domain, h) {
  
  # domain <- seq(0, 1, length.out = 100)
  
  if (domain[1] == 0) { domain[1] <- 1e-05}
  
  # 1. Generate covariance matrix 
  Sigma <- matrix(NA, length(domain), length(domain))
  h <- 0.7
  for (i in 1:length(domain)) {
    for (j in 1:length((domain))) {
      Sigma[i, j] <- frac_bm_fun(domain[i], domain[j], h)
    }
  }
  
  plot(diag(Sigma))
  # plot(Sigma[, 1])
  
  ## Simulate
  L <- t(chol(Sigma))
  z <- rnorm(length(domain), 0, 1)
  frac_bm <- L %*% z

  return(frac_bm)  
}

frac_bm_fun <- function(s, t, h) {
  0.5 * (abs(s)^(2*h) + abs(t)^(2*h) - abs(t-s)^(2*h))
}
