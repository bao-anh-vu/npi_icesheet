setwd("~/SSA_model/CNN/real_data/")

library(mvtnorm)
library(fields)
library(dplyr)
use_bed_example <- T

cond_sim_gp <- function(nsims, x_test, x_train, obs, obs_sd, use_bed_example) {

    eps <- 1e-6
    # x_test <- matrix(seq(0, 2*pi, length=n), ncol=1)
    y <- obs
    D <- rdist(x_train) 
    Dx_test <- rdist(x_test)
    Dx_test_train <- rdist(x_test, x_train)
    
    if (use_bed_example) {
        Sigma <- bed_cov(si = x_train, sj = x_train)#
        S_test_test <- bed_cov(si = x_test, sj = x_test)#
        S_test_train <- bed_cov(si = x_test, sj = x_train)
    } else {
        Sigma <- exp(-D) + diag(eps, ncol(D))
        S_test_test <- exp(-Dx_test) + diag(eps, ncol(Dx_test))
        S_test_train <- exp(-Dx_test_train) 
    }

    Si <- solve(Sigma + diag(obs_sd, ncol(Sigma))) 
    mu_p <- S_test_train %*% Si %*% y
    Sigma_p <- S_test_test - S_test_train %*% Si %*% t(S_test_train)

    sims <- rmvnorm(nsims, mu_p, Sigma_p)

    if (!is.null(obs_sd)) {
        sims_at_train <- t(rmvnorm(nsims, y, diag(obs_sd^2)))
    } else {
        sims_at_train <- matrix(rep(y, nsims), ncol=nsims)
    }

    test_df <- cbind(x_test, t(sims))
    train_df <- cbind(x_train, sims_at_train)
    sim_df <- as.data.frame(rbind(test_df, train_df))
    colnames(sim_df) <- c("x", paste0("sim", 1:nsims))

    sim_df <- dplyr::arrange(sim_df, x)

    return(list(mu_p = mu_p, Sigma_p = Sigma_p, sims = sim_df))
}

bed_cov <- function(si, sj, sill = 4000, nugget = 200, range = 50e3) { 
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

## Cov parameters
sill <- 4000
nugget <- 200
range <- 50e3

nsims <- 100
n <- 10
eps <- 1e-6

if (use_bed_example) {
    domain <- seq(0, 800e3, length.out = 100)
    x_train <- sort(sample(domain, n))
    x_test <- setdiff(domain, x_train)
    y <- 10 * sin(x_train)
    out <- cond_sim_gp(nsims, x_test, x_train, y, 
                obs_sd = rep(10, length(y)), use_bed_example = T)

} else {
    x <- seq(-0.5, 2*pi + 0.5, length=100)
    x_train <- sample(x, n) #
    # x_train <- seq(0, 2*pi, length=n)
    y <- sin(x_train)
    x_test <- setdiff(x, x_train) 
    # x_test <- seq(-0.5, 2*pi + 0.5, length=100)
    out <- cond_sim_gp(nsims, x_test, x_train, y, 0.1, use_bed_example = F)
}

# YY <- rmvnorm(100, mu_p, Sigma_p)
YY <- out$sims[, 2:(nsims+1)]
x <- out$sims[, 1]
mu_p <- out$mu_p
Sigma_p <- out$Sigma_p
q1 <- mu_p + qnorm(0.05, 0, sqrt(diag(Sigma_p)))
q2 <- mu_p + qnorm(0.95, 0, sqrt(diag(Sigma_p)))

png(paste0("./plots/bed/bed_gp.png"), width = 1000, height = 500)
matplot(x, as.matrix(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")
points(x_train, y, pch=20, cex=2)
# lines(x_test, sin(x_test), col="red")
# lines(x_test, mu_p, lwd=2)
# lines(x_test, q1, lwd=2, lty=2, col=2)
# lines(x_test, q2, lwd=2, lty=2, col=2)
dev.off()


# rbf_kernel <- function(x_test1, x_test2, length_scale = 1, variance = 1) {
#   # Computes the RBF kernel matrix between two sets of points
#   sqdist <- as.matrix(dist(rbind(x_test1, x_test2)))^2
#   K <- variance * exp(-0.5 * sqdist / length_scale^2)
#   return(K[1:nrow(x_test1), (nrow(x_test1) + 1):(nrow(x_test1) + nrow(x_test2))])
# }


# gp_regression <- function(x_test_train, y_train, x_test_test, length_scale = 1, variance = 1, noise_variance = 0.1) {
#   # Compute covariance matrices
#   K_train_train <- rbf_kernel(x_test_train, x_test_train, length_scale, variance) + noise_variance^2 * diag(nrow(x_test_train))
#   K_train_test <- rbf_kernel(x_test_train, x_test_test, length_scale, variance)
#   K_test_test <- rbf_kernel(x_test_test, x_test_test, length_scale, variance)
  
#   # Compute the inverse of K_train_train
#   K_inv <- solve(K_train_train)
  
#   # Posterior mean for test points
#   mu_test <- t(K_train_test) %*% K_inv %*% y_train
  
#   # Posterior covariance for test points
#   cov_test <- K_test_test - t(K_train_test) %*% K_inv %*% K_train_test
  
#   list(mean = mu_test, covariance = cov_test)
# }

# generate_realization <- function(mu_test, cov_test) {
#   # Draw a sample from the multivariate normal posterior
#   realization <- rmvnorm(1, as.vector(mu_test), cov_test)
#   return(realization)
# }

# # Example data
# set.seed(42)
# x_test_train <- matrix(seq(-3, 3, length.out = 10), ncol = 1)  # Training inputs
# y_train <- sin(x_test_train) + rnorm(length(x_test_train), sd = 0.1)  # Noisy training outputs

# x_test_test <- matrix(seq(-3, 3, length.out = 100), ncol = 1)   # Test inputs

# # Hyperparameters
# length_scale <- 1
# variance <- 1
# noise_variance <- 0.1

# # Run GP regression
# gp_result <- gp_regression(x_test_train, y_train, x_test_test, length_scale, variance, noise_variance)

# # Extract posterior mean and covariance for test points
# mu_test <- gp_result$mean
# cov_test <- gp_result$covariance

# # Generate a realization from the posterior
# realization <- generate_realization(mu_test, cov_test)

# # Plot the results
# plot(x_test_train, y_train, pch = 16, col = "blue", main = "GP Regression", ylim = c(-2, 2), xlab = "x_test", ylab = "f(x)")
# lines(x_test_test, sin(x_test_test), col = "green", lty = 2)  # True function
# lines(x_test_test, mu_test, col = "red")  # Posterior mean
# lines(x_test_test, realization, col = "purple", lty = 3)  # Sample from posterior

# # Optional: add uncertainty bounds (±2 standard deviations)
# std_dev <- sqrt(diag(cov_test))
# lines(x_test_test, mu_test + 2 * std_dev, col = "orange", lty = 2)
# lines(x_test_test, mu_test - 2 * std_dev, col = "orange", lty = 2)
# legend("topright", legend = c("Observed", "True Function", "Posterior Mean", "Posterior Realization", "Uncertainty Bounds"), 
#        col = c("blue", "green", "red", "purple", "orange"), lty = c(NA, 2, 1, 3, 2), pch = c(16, NA, NA, NA, NA))

