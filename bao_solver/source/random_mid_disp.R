# Random Midpoint Displacement in 1D

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

# # Example usage
# set.seed(123)
# result <- random_mid_disp(n_iterations = 11, roughness = 500)

# # Plotting
# png("./plots/random_midpoint_displacement.png", width = 800, height = 600)
# plot(result$x, result$y, type = 'l', main = "Random Midpoint Displacement", xlab = "x", ylab = "y")
# dev.off()
