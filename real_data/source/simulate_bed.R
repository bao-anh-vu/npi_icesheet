################################################################################
# This file contains functions that simulate the bedrock topography using
# conditional simulation (Cressie, 1993, page 207).
# These simulated beds are used to initialise the ensemble.
################################################################################

simulate_bed <- function(nsim, domain, obs_locations, obs, obs_sd = NULL, 
                         sill = 4000, nugget = 200, range = 50e3) {
  
  # Cov function specifications
  bed_sigma <- bed_cov(sill, nugget, range, si = domain, sj = domain)
  
  # Use the bed trend as the mean 
  # bed_mean <- create_bed(x, variable_bed = TRUE, random_bed = FALSE)
  
  df <- data.frame(obs_locations = domain[obs_locations], bed_elev = obs)
  
  # bed_fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.4,
  #                  control = loess.control(surface = "direct")) 
  bed_fit <- gam(bed_elev ~ s(obs_locations, bs = "cs", k = 20), data = df)
  bed_mean <- as.vector(predict(bed_fit, newdata = data.frame(obs_locations = domain)))
  # bed.fit <- lm(bed_elev ~ poly(obs_locations, 10), data = df) #, span = 0.12,

  ## Use loess to interpolate between observations
  ## Outside of that range just use the first and last values in the fit
  # lower_tail <- which(domain < domain[obs_locations[1]])
  # upper_tail <- which(domain > domain[obs_locations[length(obs_locations)]])
  # bed_mean <- predict(bed_fit, newdata = data.frame(obs_locations = domain))
  # bed_mean[lower_tail] <- bed_fit$fitted[1]
  # bed_mean[upper_tail] <- bed_fit$fitted[length(bed_fit$fitted)]

  # bed_mean <- bed_fit$fitted.values
  png(paste0("./plots/bed/bed_mean.png"), width = 800, height = 800)
  plot(domain, bed_mean, type = "l")
  points(domain[obs_locations], obs, col = "red")
  dev.off()

  # Unconditional simulation
  bed_sims <- matrix(0, nrow = length(domain), ncol = nsim) #matrix to store simulations
  
  L <- t(chol(bed_sigma)) 

  #  t2 <- system.time(
  for (k in 1:nsim) { #vectorise!!!
    bed_sims[, k] <- bed_mean + L %*% rnorm(length(domain))
  }
  #  )
  # bmean <- rowMeans(bed_sims)

#  png(paste0("./plots/bed/bed_mean.png"), width = 800, height = 800)
#   plot(domain, bed_sims[,1], type = "l")
#   points(domain[obs_locations], obs, col = "red")
#   dev.off()

# browser()  
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
  
## Now modify this covariance matrix so that it has the pre-specified variance at observation locations
  b.obs.sigma_mod <- b.obs.sigma
  if (!is.null(obs_sd)) {
    b.obs.sigma_mod <- b.obs.sigma + diag(obs_sd^2) # adjust the diagonal values
  #  
  } 

  sigma.inv <- solve(crossprod(chol(b.obs.sigma_mod))) #compute inverse of sigma
  
  bed.sims <- matrix(0, nrow = length(domain), ncol = nsim) # empty matrix to store simulations
  
  # t3 <- system.time(
  for (i in 1:nsim) { ## vectorise!!!
    Zns <- bed_sims[, i] 
    Zcs <- Zns + c %*% sigma.inv %*% (obs - Zns[obs_locations])
    bed.sims[, i] <- Zcs
  }

# png(paste0("./plots/bed/bed_mean.png"), width = 800, height = 800)
#   plot(domain, bed.sims[, 1], type = "l")
#   points(domain[obs_locations], obs, col = "red")
#   dev.off()

  # )
  # bed.sims.mean <- rowMeans(bed.sims) # mean of all simulated beds
  
  return(list(sims = bed.sims, mean = bed_mean))
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




