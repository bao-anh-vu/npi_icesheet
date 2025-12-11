################################################################################
# This file stores functions that create the bed and basal friction parameters.
################################################################################

## Basal friction coefficient
create_fric_coef <- function (x, L) {
  coef <- 0.02 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) 
}

## Bed topography
create_bed <- function(x, variable_bed = TRUE, random_bed = TRUE) {
  b1 <- rep(0, length(x))
  
  if (variable_bed) {
    left <- (x <= 450e3)
    right <- (x > 450e3)
    b1[left] <- - 600 + x[left] / 1000
    b1[right] <- - 150 - 5 * (x[right] / 1000 - 450)
    # b1 <- b1 + 500
  } else {
    b1 <- rep(-650, length(x))
  }
  
  ## b_r (random midpoint displacement to generate roughness of the bed)
  
  if (random_bed) {
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
    
    b1 <- b1 + br[1:length(b1)]
  }
  b1
}
