## Unconditional simulation of the basal friction coefficient

fric_cov <- function(sill, nugget, range, si, sj) {
  # matrix of distances
  d <- rdist(si, sj)
  
  # nugget model
  nug <- ifelse(d > 0, nugget, 0)
  
  # variogram
  psill <- sill - nugget
  variogram <- psill * (1 - exp(-3 * (d / range)^2)) + nug
  
  #covariance
  cov <- sill - variogram
  
  return(cov)
}

simulate_friction <- function(nsim, domain = seq(0, 800e3, by = 400), 
                              sill = 8e-5, nugget = 0, range = 2.5e3) {
  
  # Domain
  x <- domain
  
  ## Parameters
  secpera <- 31556926 #seconds per annum
  n <- 3.0 # exponent in Glen's flow law
  m <- 1/n # friction exponent
  
  ## Cov matrix for the friction
  fric.sigma <- fric_cov(sill, nugget, range, si = x, sj = x)
  
  ## Unconditional simulation of the basal friction coefficient
  Lmat <- t(chol(fric.sigma))
  fric.mu <- 0.020 #* 1e6 * (secpera)^m
  
  C_sims <- matrix(0, nrow = length(x), ncol = nsim)
  for (k in 1:nsim) { #vectorise!!!
    C_sims[, k] <- fric.mu + Lmat %*% rnorm(length(x))
  }
  
  C_sims <- C_sims #* 1e6 * secpera^m # Convert to appropriate unit (M Pa m^-1/3 s^1/3)
  
  # Cmean <- rowMeans(C_sims)
  
  return(C_sims)
}

## Simulate friction coefs
N <- 5
print("Simulating friction coefficient...")
simulated_friction <- simulate_friction(nsim = N, domain = ssa_steady$domain,
                                        sill = 8e-5, nugget = 0, range = 2500)
matplot(simulated_friction[1:350, ], type = "l") # up to GL only
