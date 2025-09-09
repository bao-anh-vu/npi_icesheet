## Compute model discrepancy

setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(mgcv)
library(dplyr)

source("./source/sim_params.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")

data_dir <- "./data/"
data_date <- "20241111" #"20241103"

## Read surface data
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

steady_state <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))

# relaxation <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_relax_", data_date, ".qs"))
vel_sim <-  relaxation$all_velocities
se_sim <- relaxation$all_top_surface
gl_sim <- relaxation$grounding_line
gl <- gl_sim[length(gl_sim)]
flowline <- relaxation$domain / 1e3
# gl_ind <- max(which(flowline <= gl))

vel_discrepancy <- vel_sim - vel_mat
se_discrepancy <- se_sim - surf_elev_mat

## Plot discrepancy
png(paste0("plots/discr/vel_discrepancy.png"), width = 800, height = 600)
matplot(flowline, vel_discrepancy, type = "l", lty = 1, col = rgb(0,0,0,0.25),
    xlab = "Distance along flowline (km)", ylab = "Velocity discrepancy (m/yr)")
abline(h = 0, col = "red", lty = 2)
dev.off()

png(paste0("plots/discr/se_discrepancy.png"), width = 800, height = 600)
matplot(flowline, se_discrepancy, type = "l", lty = 1, col = rgb(0,0,0,0.1),
    xlab = "Distance along flowline (km)", ylab = "Surface elevation discrepancy (m)")
abline(h = 0, col = "red", lty = 2)
dev.off()

## Now need to fit a Gaussian process to the discrepancy

## Assume the discrepancy in any given year is
## discrepancy = trend + fluctuations

## First average the discrepancy over time 
vel_discrepancy_avg <- rowMeans(vel_discrepancy, na.rm = TRUE)
se_discrepancy_avg <- rowMeans(se_discrepancy, na.rm = TRUE)

## Then fit a polynomial regression using lm() to the discrepancy to get the trend
vel_discr_df <- data.frame(dist = flowline, discr = vel_discrepancy_avg)
vel_discrepancy_fit <- gam(discr ~ s(dist, k = 50), data = vel_discr_df)
# vel_discrepancy_fit <- loess(vel_discrepancy_avg ~ seq(1, length(vel_discrepancy_avg)), span = 0.1)
vel_discr_trend <- predict(vel_discrepancy_fit, newdata = data.frame(dist = flowline))

se_discr_df <- data.frame(dist = flowline, discr = se_discrepancy_avg)
se_discrepancy_fit <- gam(discr ~ s(dist, k = 50), data = se_discr_df)
se_discr_trend <- predict(se_discrepancy_fit, newdata = data.frame(dist = flowline))

## Compute the fluctuations around the trend
vel_discr_fluct <- vel_discrepancy_avg - vel_discr_trend

### Plot the average discrepancy along with the polynomial fit + fluctuations
png(paste0("plots/discr/vel_discrepancy_avg.png"), width = 800, height = 600)

par(mfrow = c(2,1))
plot(flowline, vel_discrepancy_avg, type = "l",
    xlab = "Distance along flowline (km)", ylab = "Average velocity discrepancy (m/yr)")
lines(flowline, vel_discr_trend, col = "red")
abline(v = gl, lty = 2)

plot(flowline, vel_discr_fluct, type = "l", 
    xlab = "Distance along flowline (km)", ylab = "Fluctuations")
abline(h = 0, col = "red", lty = 2)
abline(v = gl, lty = 2)
dev.off()

## Now fit a Gaussian process to the fluctuations
## Might need to fit a separate GP to the grounded and floating parts of the flowline
vel_discr_fluct_df <- data.frame(dist = flowline, discr = vel_discr_fluct) %>% filter(!is.na(discr))
vel_discr_ground <- vel_discr_fluct_df %>% filter(dist <= gl)
vel_discr_float <- vel_discr_fluct_df %>% filter(dist > gl)

# RBF kernel (same helper)
rbf_kernel <- function(x1, x2, lengthscale, sigma) {
  # compute squared distances properly
#   X1 <- matrix(x1, ncol = 1)
#   X2 <- matrix(x2, ncol = 1)
  d2 <- outer(x1, x2, FUN = function(a,b) (a-b)^2)
  sigma^2 * exp(-d2 / (2 * lengthscale^2))
}

# Negative log marginal likelihood with sigma_n fixed to zero
nll_no_noise <- function(params, x, y) {
  lengthscale <- exp(params[1])
  sigma     <- exp(params[2])
  
  # Covariance (no noise) plus tiny jitter for stability
  K <- rbf_kernel(x, x, lengthscale, sigma) + 1e-8 * diag(length(x))
  
  # stable LML using Cholesky
  L <- tryCatch(t(chol(K)), error = function(e) return(NULL))
  if (is.null(L)) return(1e10) # very large penalty if K not PD
  
  # alpha <- backsolve(t(L), forwardsolve(L, y))
  v <- forwardsolve(L, y)

  log_lik <- -0.5 * t(v) %*% v - sum(log(diag(L))) - (length(x)/2) * log(2*pi)
  return(-as.numeric(log_lik))  # minimizer
}

# Example data (your vector of random numbers)
set.seed(2025)
y_g <- vel_discr_ground$discr
x_g <- vel_discr_ground$dist 

# initial guesses (log-space) for lengthscale and sigma
init <- log(c(10.0, 10.0))

opt_g <- optim(par = init,
             fn = nll_no_noise,
             x = x_g, y = y_g,
             method = "L-BFGS-B")

opt_params_g <- exp(opt_g$par)
lengthscale_g <- opt_params_g[1]
sigma_g     <- opt_params_g[2]

cat("Optimized (noise-free) hyperparameters:\n")
cat("Lengthscale (ground) =", lengthscale_g, "\n")
cat("SD (ground) =", sigma_g, "\n")

y_f <- vel_discr_float$discr #%>% filter(!is.na(discr))
x_f <- vel_discr_float$dist

# initial guesses (log-space) for lengthscale and sigma
init <- log(c(10.0, 10.0))
# test <- nll_no_noise(init, x_f, y_f)  # check it runs
# print(test)
# browser()

opt_f <- optim(par = init,
             fn = nll_no_noise,
             x = x_f, y = y_f,
             method = "L-BFGS-B")

opt_params_f <- exp(opt_f$par)
lengthscale_f <- opt_params_f[1]
sigma_f   <- opt_params_f[2]

cat("Optimized (noise-free) hyperparameters:\n")
cat("Lengthscale (float)  =", lengthscale_f, "\n")
cat("Sd (float) =", sigma_f, "\n")


browser()

# Posterior prediction (noise-free)
x_new <- seq(1, length(y), length.out = 200)
K    <- rbf_kernel(x, x, lengthscale, sigma) + 1e-8 * diag(length(x))
K_s  <- rbf_kernel(x, x_new, lengthscale, sigma)
K_ss <- rbf_kernel(x_new, x_new, lengthscale, sigma) + 1e-8 * diag(length(x_new))

# Solve using Cholesky for stability
L <- chol(K)
K_inv_y <- backsolve(t(L), forwardsolve(L, y))

mu_s  <- t(K_s) %*% K_inv_y
# predictive covariance
V <- forwardsolve(L, K_s)
cov_s <- K_ss - t(V) %*% V
std_s <- sqrt(pmax(0, diag(cov_s)))  # clamp small negatives due to numerical error

# Plot: GP will interpolate training points exactly
plot(x, y, pch = 16, ylim = range(c(y, mu_s + 2*std_s, mu_s - 2*std_s)),
     main = "GP (noise-free) — interpolating fit", xlab = "x", ylab = "y")
lines(x_new, mu_s, lwd = 2)
polygon(c(x_new, rev(x_new)),
        c(mu_s + 2*std_s, rev(mu_s - 2*std_s)),
        col = rgb(0,0,1,0.15), border = NA)

# Optionally check posterior variance at training points (should be ~0)
train_var <- diag(cov_s)[match(x, round(x_new,6))]  # approximate indexing
print(head(train_var))

