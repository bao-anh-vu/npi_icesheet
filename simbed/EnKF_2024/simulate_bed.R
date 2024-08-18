################################################################################
# This file contains functions that simulate the bedrock topography using
# conditional simulation (Cressie, 1993, page 207).
# These simulated beds are used to initialise the ensemble.
################################################################################

simulate_bed <- function(nsim, domain, obs_locations, obs, 
                         sill = 4000, nugget = 200, range = 50e3) {
  # Cov function specifications
  bed_sigma <- bed_cov(sill, nugget, range, si = domain, sj = domain)
  
  # Use the bed trend as the mean 
  # bed_mean <- create_bed(x, variable_bed = TRUE, random_bed = FALSE)
  
  df <- data.frame(obs_locations = domain[obs_locations], bed_elev = obs)
  bed.fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.25, 
                   control = loess.control(surface = "direct")) 
  bed_mean <- predict(bed.fit, newdata = data.frame(obs_locations = domain))
  
  # Unconditional simulation
  bed_sims <- matrix(0, nrow = length(domain), ncol = nsim) #matrix to store simulations
  
  for (k in 1:nsim) { #vectorise!!!
    bed_sims[, k] <- bed_mean + t(chol(bed_sigma)) %*% rnorm(length(domain))
  }
  bmean <- rowMeans(bed_sims)
  
  ## Plot unconditional bed simulations
  # yrange <- c(-1500, max(bed_sims))
  # plot(x[500:1500] / 1000, bed_sims[500:1500, 1], type = "l", ylim = yrange, col = "lavender",
  #      xlab = "x (km)", ylab = "Bed elevation")
  # lines(x[500:1500] / 1000, bed_sims[500:1500, 2], type = "l", col = "lightblue")
  # lines(x[500:1500] / 1000, bmean[500:1500], type = "l", col = "black")
  # 
  # legend("bottomleft", legend = c("1st member", "2nd member", "Ensemble mean", "True bed"),
  #        lty = c(1, 1, 1, 2), col = c("lavender", "lightblue", "black", "red"), cex = 0.5)
  
  ## Conditional simulation (see Cressie 1993 pg. 207 for details)
  
  # Covariance between all grid points and observation locations
  c <- bed_cov(sill, nugget, range, si = domain, sj = domain[obs_locations])
  
  # Covariance between observation locations
  b.obs.sigma <- bed_cov(sill, nugget, range, domain[obs_locations], domain[obs_locations])
  
  sigma.inv <- solve(crossprod(chol(b.obs.sigma))) #compute inverse of sigma
  
  bed.sims <- matrix(0, nrow = length(domain), ncol = nsim) # empty matrix to store simulations
  
  for (i in 1:nsim) {
    Zns <- bed_sims[, i] 
    Zcs <- Zns + c %*% sigma.inv %*% (obs - Zns[obs_locations])
    bed.sims[, i] <- Zcs
  }
  
  bed.sims.mean <- rowMeans(bed.sims) # mean of all simulated beds
  
  # par(mfrow = c(1,1))
  # # plot(x/1000, b, type = "l")
  # plot(x/1000, bed.sims[, 1], type = "l", col = "lightblue") # first member
  # lines(x/1000, bed.sims.mean, lty = 2, col = "red") # conditional sim
  # points(x[obs_locations]/1000, obs, col = "cyan", bg = "cyan", pch = 20)
  # legend("topright", legend = c("True bed", "Mean simulated bed", "First ensemble member"),
  #        col = c("black", "red", "lightblue"), lty = 1:2, cex = 0.5)
  # 
  return(bed.sims)
}

bed_cov <- function(sill, nugget, range, si, sj) { 
  # matrix of distances
  d <- rdist(si, sj)

  # nugget model
  nug <- ifelse(d > 0, nugget, 0)

  # variogram
  psill <- sill - nugget
  variogram <- psill * (1 - exp(-3 * d / range)) + nug

  #covariance
  cov <- sill - variogram
  # s * (1 - exp(-3 * spDists(coordinates(x),coordinates(y)) / range))
  return(cov)
}


# plot(bed_sigma[1, ])




