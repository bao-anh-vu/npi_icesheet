## Test velocities ##
library("parallel")
library("Matrix")
library("qlcMatrix")
library("R.utils")
library("mvtnorm")
library("fields")
library("dplyr")
library("matrixStats")

source("surface_elev.R") # for calculating grounding line position
source("solve_thickness.R") # solving for the ice thickness
source("solve_velocity_azm.R") # solving for the velocity
source("azm_cond_sim.R") # for generating correlated noise
source("propagate.R") # for propagating ensemble
source("obs_operator.R") # observation operator 

## Flags
add_process_noise <- TRUE # add process noise every year
use_cov_taper <- TRUE # use covariance tapering

## Ensemble parameters
N <- 20 # number of ensemble members
years <- 1
steps_per_yr <- 100

##### I. Setup #####
print("Loading data...")

# Reference ice thickness, velocity, surface elevation etc.
reference <- readRDS(file = "/home/babv971/SSA_model/EnKF/Output/reference_20220202.rds")

# Observations of surface velocity and surface elevation
observations <- readRDS(file = "/home/babv971/SSA_model/EnKF/Output/observations_20220202.rds")

domain <- reference$domain 
J <- length(domain) # number of grid points

# Initial ensemble
ini_ens <- readRDS(file = "/home/babv971/SSA_model/EnKF/Output/ini_ens_20220124.rds")

# Initial velocity
ini_velocity <- matrix(rep(as.vector(reference$ini_velocity), N), J, N) 
mean_prev_velocity <- rowMeans(ini_velocity)
prev_velocity <- ini_velocity

## Set up ensemble
ens <- ini_ens
prev_velocity <- ini_velocity

## Covariance tapering setup
wendland <- function(theta, D) { ## Wendland function
  R <- D / theta
  W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
}

taper_mat <- wendland(0.1 * domain[length(domain)], rdist(domain, domain))
ens_taper1 <- kronecker(matrix(rep(1, 6), 3, 2), taper_mat) # for the PH term
ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat) # for the HPH term

##### II. Propagation #####
print("Propagating ensemble...")
mc.t1 <- proc.time()

# Append the velocities to the ensemble
ens_all <- rbind(ini_ens, prev_velocity)

# Convert ensemble from matrix to list
ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])


# Apply the propagate function to every ensemble member
ens.list <- mclapply(ens.list, propagate, mc.cores = 6L, 
                     domain = domain, steps_per_yr = steps_per_yr)

# Convert ensemble from list back to matrix
ens_all <- matrix(unlist(ens.list), nrow = 4*J, ncol = N)

## Extract forecast ens and velocity ens
forecast_ens <- ens_all[1:(3*J), ]
prev_velocity <- ens_all[(3*J+1):(4*J), ]

# Save the mean velocity
# mean_prev_velocity <- rowMeans(prev_velocity)

if (add_process_noise) { ## Add process noise
  print("Adding process noise...")
  h_sd <- pmin(0.02 * rowMeans(forecast_ens[1:J, ]), 20)
  
  h_noise <- cond_sim(mu = rep(0, J), sd = h_sd,
                      domain, l = 50e3, nsims = N)
  
  forecast_ens[1:J, ] <- as.matrix(forecast_ens[1:J, ] + h_noise[, 1:N]) # add noise to ice thickness
  
}

##### III. Analysis #####
print("Running analysis...")

analysis.t1 <- proc.time()

# Compute P %*% t(H) and H %*% P %*% t(H)
# HPH <- NULL #matrix(NA, nrow = 3*J, ncol = 3*J)
# PH <- NULL #matrix(NA, nrow = 3*J, ncol = 2*J)

## Apply observation operator to every ensemble member
HX <- mclapply(ens.list, obs_operator,
               domain = domain, mc.cores = 6L)

HX <- matrix(unlist(HX), nrow = 2*J, ncol = N)

# Compute P %*% t(H) and H %*% P %*% t(H)
HPH <- 1 / (N - 1) * tcrossprod(HX - rowMeans(HX))
PH <- 1 / (N - 1) * tcrossprod(forecast_ens - rowMeans(forecast_ens), 
                               HX - rowMeans(HX))

## Construct measurement error covariance matrix R (depends on the velocity)
vel_std <- pmin(0.25 * rowMeans(HX)[(J+1):(2*J)], 20)
vel_std[vel_std == 0] <- 1e-05
R <- diag(c(rep(10^2, J), vel_std^2)) # measurement noise for top surface and velocity

## Calculate Kalman gain matrix K 
print("Calculating Kalman gain matrix...")
if (use_cov_taper) {
  
  ##### CHANGED HERE #######
  # Tapered Kalman gain
  Sigma <- as(ens_taper2 * HPH + R, "dsCMatrix")
  #L <- t(chol(ens_taper2 * HPH + R))
  #Linv <- solve(L)
  #K <- (ens_taper1 * PH) %*% t(Linv) %*% Linv #solve(hadamard(ens_taper2, HPH) + R)
  
} else {
  
  # K_test <- PH %*% solve(HPH + R) # compute the original Kalman gain matrix for comparison
  L <- t(chol(HPH + R))
  K <- PH %*% solve(t(L)) %*% solve(L)
  
}

print("Updating ensemble...")
## Extract yearly observations
year_ind <- 2 # t = 1a
y <- c(observations$surface_obs[, year_ind], observations$velocity_obs[, year_ind])

## Generate "observation ensemble" Y = (y, y, ..., y)
Y <- matrix(rep(y, N), 2*J, N) 

## Update ensemble

##### CHANGED HERE #######
#noise <- rmvnorm(N, rep(0, 2*J), R)
noise <- matrix(rnorm(N*2*J, mean = 0, sd = sqrt(diag(R))), nrow = N)

##### CHANGED HERE #######
#ens <- forecast_ens + K %*% (Y - HX - t(noise))
ens <- forecast_ens + (ens_taper1 * PH) %*% solve(a = Sigma, b = (Y - HX - t(noise)))

# Store analysis ensemble and ensemble mean
# enkf_means[, year_ind] <- rowMeans(ens)
# enkf_ens[[year_ind]] <- ens

analysis.t2 <- proc.time()

##### IV. Update velocity at the 100th step based on the estimated thickness #####
ens_all <- rbind(ens, prev_velocity) # append velocity to ensemble

# Convert ensemble to list
ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])

# Propagate velocity
velocity.list <- mclapply(ens.list, 
                          function(v) { 
                            # evolve the velocity
                            u <- as.vector(solve_velocity(prev_velocity = v[(3*J+1):(4*J)], #mean_prev_velocity,
                                                          thickness = v[1:J], 
                                                          domain = domain,
                                                          bed = v[(J+1):(2*J)],
                                                          friction = 10^v[(2*J+1):(3*J)],
                                                          perturb_hardness = TRUE))
                            
                            return(u)
                            
                          }, 
                          mc.cores = 6L)

# Convert velocity list to matrix
prev_velocity <- matrix(unlist(velocity.list), nrow = J, ncol = N)

mc.t2 <- proc.time()

################################ Plots #########################################
par(mfrow = c(1,1))

plot(reference$all_thicknesses[, 2], type = "l") # True ice thickness at t = 1a
lines(rowMeans(ens[1:J, ]), col = "red") # Filtered ice thickness

plot(reference$all_velocities[, 2], type = "l") # True velocity at t = 1a
lines(rowMeans(prev_velocity), col = "red")

## Benchmarking
print(mc.t2 - mc.t1)
