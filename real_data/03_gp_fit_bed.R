## Fit a Gaussian process to bed elevation data along the flowline
## to use as a prior

setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(dplyr)
library(mvtnorm)
library(fields)

## Flags
cov_fun <- "exp" # or "rbf"
use_bed_data <- T

if (use_bed_data) {
  ## Flowline
  flowline <- qread(paste0("./data/flowline_regrid.qs"))
  J <- nrow(flowline)
  # J <- 2001 # number of grid points
  # flowline <- flowline[1:J, ]
  flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J - 1)])^2 + (flowline$y[2:J] - flowline$y[1:(J - 1)])^2)
  flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

  ## Bed data
  bed <- qread(file = paste0("./data/bedmap_obs.qs"))
  bed_sd <- unlist(qread(file = paste0("./data/bedmap_sd.qs")))
  bed <- bed[1:J]
  bed_sd <- bed_sd[1:J]

  bed_df <- data.frame(ind = 1:J, loc = flowline_dist, bed_elev = bed, bed_sd = bed_sd)
  bed_obs_df <- bed_df %>% filter(!is.na(bed_elev) & !is.na(bed_sd)) # %>% pull(ind)

qsave(bed_obs_df, file = "./data/bedmap/bed_obs_df_all.qs")

  x <- bed_obs_df$loc
  y <- bed_obs_df$bed_elev
  xgrid <- bed_df$loc # seq(0, 1, length.out = 200)

## Initial params for optim
sigma <- as.numeric(bed_obs_df$bed_sd)
  if (cov_fun == "exp") {
    init_params <- log(c(s = 1, ls = 20e3)) #, sigma = 10))
  } else {
    init_params <- log(c(s = 100, ls = 10e3)) #, sigma = 10))
  }
} else {
  # Training data
  # x <- c(0, 0.25, 0.5, 0.75, 0.8, 1)
  # y <- c(0.1, 0.5, 0.9, 0.6, 0.5, 0.2)

  # set.seed(2025)
  x <- seq(0, 1, length.out = 10)
  y <- rnorm(length(x), sd = 1)

  # Prediction grid
  xgrid <- seq(0, 1, length.out = 100)

  ## Initial params for optim
  sigma <- seq(0.01, 0.05, length.out = length(x)) # initial noise level
  if (cov_fun == "exp") {
    init_params <- log(c(s = 1, ls = 0.1)) #, sigma = 0.1))
  } else {
    init_params <- log(c(s = 1, ls = 0.1)) #, sigma = 0.1))
  }
  
  
}

# Design matrix: constant mean -> just a column of ones
n <- length(x)
H <- matrix(1, n, 1)

# --- Exp covariance function ---
exp_kernel <- function(si, sj, s, lengthscale) {
  # matrix of distances
  d <- rdist(si, sj)
  cov <- s^2 * exp(-d / lengthscale)
  return(cov)
}

# --- RBF kernel ---
rbf_kernel <- function(x1, x2, s, lengthscale) {
  # X1 <- matrix(x1, ncol = 1)
  # X2 <- matrix(x2, ncol = 1)
  # d2 <- as.matrix(dist(rbind(X1, X2)))^2
  # n1 <- nrow(X1)
  # n2 <- nrow(X2)
  # d2_block <- d2[1:n1, (n1 + 1):(n1 + n2)]
  # cov <- s^2 * exp(-0.5 * d2_block / (lengthscale^2))
  
  d <- rdist(x1, x2)
  cov <- s^2 * exp(-0.5 * (d / lengthscale)^2)

  return(cov)
}

# --- Stable linear algebra helpers ---
chol_solve <- function(Sigma, b) {
  R <- chol(Sigma)
  w <- forwardsolve(t(R), b)
  backsolve(R, w)
}
# chol_logdet <- function(Sigma) {
#   R <- chol(Sigma)
#   2 * sum(log(diag(R)))
# }

# --- Negative log marginal likelihood ---
neg_log_marginal <- function(params, sigma, cov_fun, X, y, H) {
  s <- exp(params[1])
  ls <- exp(params[2])
  # sigma <- exp(params[3])
  n <- length(y)

  if (cov_fun == "exp") {
    K <- exp_kernel(X, X, s, ls)
  } else {
    K <- rbf_kernel(X, X, s, ls)
  }

  Sigma <- K + (sigma^2) * diag(n)

  Sinv_y <- chol_solve(Sigma, y) # compute Sigma^{-1} y
  Sinv_H <- chol_solve(Sigma, H) # compute Sigma^{-1} H

  A <- t(H) %*% Sinv_H
  bvec <- t(H) %*% Sinv_y
  beta_hat <- solve(A, bvec)

  res <- y - H %*% beta_hat
  quad <- t(res) %*% chol_solve(Sigma, res) # compute (y - H*beta)^T Sigma^{-1} (y - H*beta)
  # logdet <- chol_logdet(Sigma)
  R <- chol(Sigma)  
  

  nll <- 0.5 * (quad + 2 * sum(log(diag(R))) + n * log(2 * pi))
  cat("nll = ", as.numeric(nll), " for params: ", round(c(s, ls), 4), "\n")
  
  return(as.numeric(nll))
}

# --- Fit by optimizing hyperparameters ---

opt <- optim(
  par = init_params,
  fn = neg_log_marginal,
  X = x, y = y, H = H, sigma = sigma, cov_fun = cov_fun,
  method = "L-BFGS-B",
  lower = log(c(1e-6, 1e-6, 1e-6)),
  upper = log(c(1e6, 1e6, 1e6))
)

fitted <- exp(opt$par)
names(fitted) <- c("s", "lengthscale")#, "sigma_noise")
print(fitted)

# --- Compute beta_hat at optimum ---
s <- fitted["s"]
ls <- fitted["lengthscale"]
# sigma <- fitted["sigma_noise"]

if (cov_fun == "exp") {
  K <- exp_kernel(x, x, s, ls)
} else {
  K <- rbf_kernel(x, x, s, ls)
}

Sigma <- K + sigma^2 * diag(n)

Sinv_y <- chol_solve(Sigma, y)
Sinv_H <- chol_solve(Sigma, H)
A <- t(H) %*% Sinv_H
beta_hat <- solve(A, t(H) %*% Sinv_y) # (H^T Sigma^{-1} H)^{-1} H^T Sigma^{-1} y
cat("Estimated mean (beta0):", beta_hat, "\n")

# A_inv <- solve(A)

# --- Prediction function ---
predict_gp <- function(xstar, x, cov_fun = "exp", beta_hat, s, ls, sigma) {
  n <- length(x)
  nstar <- length(xstar)
  Hstar <- matrix(1, nstar, 1) # constant mean
  # kstar <- t(sapply(xstar, function(xs) exp_kernel(xs, x, s, ls)))
  # Kss <- outer(xstar, xstar, function(a,b) exp_kernel(a, b, s, ls))
  if (cov_fun == "exp") {
    K <- exp_kernel(x, x, s, ls)
  } else {
    K <- rbf_kernel(x, x, s, ls)
  }
  Sigma <- K + sigma^2 * diag(n)
  Sigma_inv <- solve(Sigma)

  if (cov_fun == "exp") {
    kstar <- exp_kernel(xstar, x, s, ls)
    Kss <- exp_kernel(xstar, xstar, s, ls)
  } else {
    kstar <- rbf_kernel(xstar, x, s, ls)
    Kss <- rbf_kernel(xstar, xstar, s, ls)
  }

  # Sinv <- function(v) chol_solve(Sigma, v)
  res <- y - H %*% beta_hat
  mu_star <- Hstar %*% beta_hat + kstar %*% chol_solve(Sigma, res)
  cov_star <- Kss - kstar %*% solve(Sigma) %*% t(kstar)

  # mu_star <- Hstar %*% beta_hat + kstar %*% Sinv(res)

  # v_obj <- t(apply(kstar, 1, function(krow) chol_solve(Sigma, krow)))
  # term2 <- rowSums(kstar * v_obj)
  # cov_star <- Kss - term2


  # # add correction for beta uncertainty
  # Ht_Sinv_k <- t(H) %*% t(v_obj)

  # diff <- t(Hstar) - Ht_Sinv_k
  # cov_correction <- t(diff) %*% A_inv %*% diff
  # cov_star <- cov_star + cov_correction

  list(mean = mu_star, cov = cov_star)
}

# --- Test predictions and plot ---

pred <- predict_gp(
  cov_fun = cov_fun,
  xstar = xgrid, x = x,
  beta_hat = beta_hat,
  s = fitted["s"], ls = fitted["lengthscale"], sigma = sigma#fitted["sigma_noise"]
)
upper <- pred$mean + 2 * sqrt(pmax(0, diag(pred$cov)))
lower <- pred$mean - 2 * sqrt(pmax(0, diag(pred$cov)))

L <- t(chol(pred$cov))
sim1 <- pred$mean + L %*% rnorm(nrow(L))
sim2 <- pred$mean + L %*% rnorm(nrow(L))

if (use_bed_data) {
  
  plot_name <- "GP_fit_to_bed.png"
} else {
  plot_name <- "GP_fit_to_toy.png"
}

png(paste0("./plots/bed/", plot_name), width = 800, height = 500)
plot(xgrid, pred$mean,
  type = "l", lwd = 2, col = "blue",
  ylim = range(c(upper, lower)), 
  xlab = "x", ylab = "y", main = "GP fit" 
)
lines(xgrid, sim1, lwd = 2, col = "goldenrod")
lines(xgrid, sim2, lwd = 2, col = "salmon")
if (use_bed_data) {
  lines(xgrid, bedmachine$bed_avg, lwd = 2, col = "forestgreen")
}
lines(xgrid, upper, lty = 2, col = "red")
lines(xgrid, lower, lty = 2, col = "red")
points(x, y, pch = 19, col = "black")
legend("topright",
  legend = c("data", "pred mean", "±2 sd"),
  pch = c(19, NA, NA), lty = c(NA, 1, 2), col = c("black", "blue", "red")
)
dev.off()

GP_fit <- list(beta = beta_hat, cov_fun = cov_fun, var = s, lengthscale = ls, 
              mean = pred$mean, cov = pred$cov, x = x, y = y, xstar = xgrid)

# qsave(fitted, file = paste0("./data/bed/GP_fit_params_", cov_fun, ".qs"))
qsave(GP_fit, file = paste0("./data/bedmap/GP_fit_", cov_fun, ".qs"))
