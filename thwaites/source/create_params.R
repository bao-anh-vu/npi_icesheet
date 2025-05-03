################################################################################
# This file stores functions that create the bed and basal friction parameters.
################################################################################

## Basal friction coefficient
create_fric_coef <- function (x, L) {
  coef <- 0.02 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) 
}

## Bed topography
create_bed <- function(x, variable_bed = TRUE, random_bed = TRUE) {
  
  b_trend <- ifelse(x <= 450e3, -600 + x / 1e3, -150 - 5 * (x / 1e3 - 450))
  b_noise <- random_mid_disp(n_iterations = 12, roughness = 500)

  bed <- b_trend + b_noise$y[1:length(x)]

  return(bed)

  # ------------ Old code ------------
  # b1 <- rep(0, length(x))
  
  # if (variable_bed) {
  #   left <- (x <= 450e3)
  #   right <- (x > 450e3)
  #   b1[left] <- - 600 + x[left] / 1000
  #   b1[right] <- - 150 - 5 * (x[right] / 1000 - 450)
  #   # b1 <- b1 + 500
  # } else {
  #   b1 <- rep(-650, length(x))
  # }
  
  ## b_r (random midpoint displacement to generate roughness of the bed)
  
  # if (random_bed) {
  #   K <- 2
  #   br <- rep(0, K)
  #   mu <- 0
  #   sigma <- 500
  #   h <- 0.7
  #   reps <- 11
    
  #   for (i in 1:reps) {
  #     midpoints <- 0.5 * (br[2:length(br)] + br[1:(length(br) - 1)])
  #     midpoints <- midpoints + rnorm(length(midpoints), mu, sigma)
  #     indices <- seq(2, length(br), 1)
  #     br <- insert(br, indices, midpoints)
  #     sigma <- sigma / (2^h)
  #   }
    
  #   b1 <- b1 + br[1:length(b1)]
  # }
  # b1
}


random_mid_disp <- function(n_iterations, roughness, seed = 123) {
  x <- c(0, 1) # Start with two points at 0 and 1
  y <- c(0, 0) # y-coordinates of the starting points
  
    set.seed(seed) # Set seed for reproducibility

  for (i in 1:n_iterations) {
    # Insert midpoints between each pair of existing points
    mid_x <- (x[-length(x)] + x[-1]) / 2
    mid_y <- (y[-length(y)] + y[-1]) / 2 + rnorm(length(mid_x), mean = 0, sd = roughness)
    
    # Merge x and mid_x
    x <- as.vector(rbind(x[-length(x)], mid_x))
    x <- c(x, 1)
    
    # Merge y and mid_y
    y <- as.vector(rbind(y[-length(y)], mid_y))
    y <- c(y, 0)
    
    roughness <- roughness / 2#^0.7  # Decrease roughness
  }
  
  list(x = x, y = y)
}