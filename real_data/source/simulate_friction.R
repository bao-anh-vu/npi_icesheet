## Unconditional simulation of the basal friction coefficient

# library("gstat")
# library("fields")
# library("sp")
# library("R.utils")
# 
# rm(list = ls())
# set.seed(123)

# ## Parameters
# secpera <- 31556926 #seconds per annum
# n <- 3.0 # exponent in Glen's flow law
# m <- 1/n # friction exponent
# 
# ## Ice shelf specifications
# L = 800e3 # length of domain
# 
# # Discretise domain
# x0 <- 0
# J <- 2000 # number of steps
# dx <- L / J # increments
# x <- seq(x0, L, dx)
# 
# ## Basal friction coefficient
# fric_coef <- function (x, L) {
#   coef <- 0.025 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) # Basal friction coefficient
#   coef
# }
# C <- fric_coef(x, L) * 1e6 * (secpera)^m
# 
# plot(x[1:500], C[1:500], type = "l", xlab = "x (m)",
#      ylab = "Basal friction coef (M Pa m-1/3 s1/3)")

### Simulation ###


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

fric_cov2 <- function(sigma, scale, si, sj) { # use a covariance function directly without variogram
  # squared exponential function
  d <- rdist(si, sj)
  cov <- sigma^2 * exp(-3 * (d/scale)^2)
  diag(cov) <- diag(cov) + 1e-8
  
  return(cov)
} 

simulate_friction2 <- function(nsim, domain, sill = 8e-5, nugget = 0, range = 2.5e3) {
  
  # Domain
  x <- domain
  
  # Covariance function specifications
  # sill <- 8e-5 # pass these into the function as arguments
  # nugget <- 0
  # range <- 25e3
  
  ## Parameters
  secpera <- 31556926 #seconds per annum
  n <- 3.0 # exponent in Glen's flow law
  m <- 1/n # friction exponent
  
  # xy <- expand.grid(domain, domain)
  # names(xy) <- c("x","y")
  # gridded(xy) = ~x+y
  
  # test <- gstat::vgm(model = "Exp", psill = sill - nugget, range = range)
  # g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,
  #                  model = test, nmax = 10)
  # yy <- predict(g.dummy, xy, nsim = 2)
  
  ## Cov matrix for the friction
  # fric.sigma <- fric_cov(sill, nugget, range, si = x, sj = x)
  fric.sigma <- fric_cov2(sqrt(sill), range, si = x, sj = x)
  
  ## Unconditional simulation of the basal friction coefficient
  Lmat <- t(chol(fric.sigma))
  fric.mu <- 0.025 #* 1e6 * (secpera)^m
  
  C_sims <- matrix(0, nrow = length(x), ncol = nsim)
  for (k in 1:nsim) { #vectorise!!!
    C_sim <- fric.mu + Lmat %*% rnorm(length(x))
    
    if (any(C_sim < 0)) {
      C_sim <- ifelse(C_sim > 0, C_sim, 0.001) # in case friction value is negative, replace with a small value
    }
    
    C_sims[, k] <- C_sim
  }
  
  C_sims <- C_sims * 1e6 * secpera^m # Convert to appropriate unit (M Pa m^-1/3 s^1/3)
  
  Cmean <- rowMeans(C_sims)
  
  # ### Plot simulations ##
  # return(list(fric_cov = fric.sigma, fric_sims = C_sims))
  return(C_sims)
  
}

simulate_friction <- function(nsim, domain, sill = 8e-5, nugget = 0, range = 2.5e3) {
  
  # Domain
  x <- domain
  
  # Covariance function specifications
  # sill <- 8e-5 # pass these into the function as arguments
  # nugget <- 0
  # range <- 25e3
  
  ## Parameters
  secpera <- 31556926 #seconds per annum
  n <- 3.0 # exponent in Glen's flow law
  m <- 1/n # friction exponent
  
  # xy <- expand.grid(domain, domain)
  # names(xy) <- c("x","y")
  # gridded(xy) = ~x+y
  
  # test <- gstat::vgm(model = "Exp", psill = sill - nugget, range = range)
  # g.dummy <- gstat(formula = z~1, dummy = TRUE, beta = 0,
  #                  model = test, nmax = 10)
  # yy <- predict(g.dummy, xy, nsim = 2)
  
  ## Cov matrix for the friction
  fric.sigma <- fric_cov(sill, nugget, range, si = x, sj = x)
  # fric.sigma <- fric_cov2(sqrt(sill), range, si = x, sj = x)
  
  ## Unconditional simulation of the basal friction coefficient
  Lmat <- t(chol(fric.sigma))
  fric.mu <- 0.020 #* 1e6 * (secpera)^m
  
  C_sims <- matrix(0, nrow = length(x), ncol = nsim)
  for (k in 1:nsim) { #vectorise!!!
    C_sim <- fric.mu + Lmat %*% rnorm(length(x))
    
    if (any(C_sims < 0)) {
      C_sim[C_sim < 0] <- 0.001
    }
    
    C_sims[, k] <- C_sim
  }
  
  C_sims <- C_sims * 1e6 * secpera^m # Convert to appropriate unit (M Pa m^-1/3 s^1/3)
  
  Cmean <- rowMeans(C_sims)
  
  # ### Plot simulations ##
  return(C_sims)
}
