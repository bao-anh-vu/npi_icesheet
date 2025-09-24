################################################################################
# This file stores functions that create the bed and basal friction parameters.
################################################################################

## Basal friction coefficient
create_fric_coef <- function (x, L) {
  coef <- 0.02 + 0.015 * sin(1 * 2 * pi * x / L) * sin(32 * pi * x / L)
  return(coef)
}

## Bed topography
create_bed <- function(x, bed_mean = NULL, random_bed = TRUE) {
  # bed_mean <- rep(0, length(x))
  
  if (!is.null(bed_mean)) {
    bed_mean <- bed_mean
  } else {
    left <- (x <= 150e3)
    right <- (x > 150e3)
    bed_mean[left] <- - 1000 + x[left] / 1000
    bed_mean[right] <- - 850 - 2 * (x[right] / 1000 - 150)
    # bed_mean <- bed_mean - 500
  } #else {
  #   bed_mean <- rep(-650, length(x))
  # }

  ## b_r (random midpoint displacement to generate roughness of the bed)
  
  if (random_bed) {
    set.seed(2024)
    K <- 2
    br <- rep(0, K)
    mu <- 0
    sigma <- 500
    h <- 0.7
    reps <- 11
    
    for (i in 1:reps) {
      midpoints <- 0.5 * (br[2:length(br)] + br[1:(length(br) - 1)])
      midpoints <- midpoints + rnorm(length(midpoints), mu, sigma)
      indices <- seq(2, length(br), 1)
      br <- insert(br, indices, midpoints)
      sigma <- sigma / (2^h)
    }
    
    bed_mean <- bed_mean + br[1:length(bed_mean)]
  }
  bed_mean
}

