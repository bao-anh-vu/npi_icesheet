## Posterior predictive distribution for real data examplesetwd("~/SSA_model/CNN/real_data/")

setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
library(R.utils)
# library("sp")
library(fields)
# library("tidyr")
# library(dplyr)
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(abind)
# library(ggplot2)
# library(gridExtra)
# library(FRK)
library(qs)
library(parallel)

source("./source/sim_obs.R")
# source("./source/cond_sim_gp.R")
source("./source/get_ini_thickness.R")
source("./source/process_sim_results.R")
# source("./source/fit_basis.R")
source("./source/surface_elev.R")
# source("./source/create_params.R")
# source("./source/create_ref.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/get_surface_obs.R")
# source("./source/simulate_bed.R")
# source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

data_date <- "20241111" # "20241103"
sets <- 51:100 # 6:20
use_missing_pattern <- T
use_basal_melt_data <- T
correct_model_discrepancy <- T

setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

if (use_missing_pattern) {
  data_dir <- paste0("./data/training_data/", setsf, "/missing/")
  output_dir <- paste0("./output/cnn/", setsf, "/missing/")
  plot_dir <- paste0("./plots/cnn/", setsf, "/missing/")
} else {
  data_dir <- paste0("./data/training_data/", setsf, "/nonmissing/")
  output_dir <- paste0("./output/cnn/", setsf, "/nonmissing/")
  plot_dir <- paste0("./plots/cnn/", setsf, "/nonmissing/")
}

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")


## Posterior samples
post_fric_samples <- qread(file = paste0(output_dir, "fric_samples_real_", data_date, ".qs"))
post_bed_samples <- qread(file = paste0(output_dir, "bed_samples_real_", data_date, ".qs"))

fric_post_mean <- qread(file = paste0(output_dir, "pred_fric_real_", data_date, ".qs"))
bed_post_mean <- qread(file = paste0(output_dir, "pred_bed_real_", data_date, ".qs"))


## Compare with simulations from the prior
set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")
prior_fric_samples <- t(qread(file = paste0("./data/training_data/friction_arr_", setf, "_", data_date, ".qs")))
prior_bed_samples <- t(qread(file = paste0("./data/training_data/bed_arr_", setf, "_", data_date, ".qs")))

png(file = paste0(plot_dir, "/posterior_vs_prior_", data_date, ".png"), width = 800, height = 800)
matplot(post_fric_samples[, 1:3], type = "l", col = "red", main = "Friction samples vs prior simulations", ylab = "Friction coeff (scaled)", xlab = "Sample index")
matlines(prior_fric_samples[, 1:3], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)

matplot(post_bed_samples[, 1:3], type = "l", col = "red", main = "Bed samples vs prior simulations", ylab = "Bedrock (m)", xlab = "Sample index")
matlines(prior_bed_samples[, 1:3], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)
dev.off()

n_post_samples <- 100
post_param_list <- lapply(1:n_post_samples, function(r) {
  list(
    friction = post_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = post_bed_samples[, r] #+ bed_mean
  )
})

prior_param_list <- lapply(1:n_post_samples, function(r) {
  list(
    friction = prior_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = prior_bed_samples[, r] #+ bed_mean
  )
})

# param_mean <- list()
# param_mean[[1]] <- list(
#     friction = as.vector(fric_post_mean), #* 1e6 * params$secpera^(1 / params$n),
#     bedrock = as.vector(bed_post_mean) #+ bed_mean
# )

## Plot parameters
png(filename = paste0(plot_dir, "fric_bed_samples_real_", data_date, ".png"), width = 800, height = 800)
par(mfrow = c(2, 1))
plot(post_param_list[[1]]$friction, type = "l", main = "Friction sample", ylab = "Friction coeff (scaled)", xlab = "Sample index")
plot(post_param_list[[1]]$bedrock, type = "l", main = "Bed sample", ylab = "Friction coeff (scaled)", xlab = "Sample index")
dev.off()

############################################################
##      Simulate observations based on posterior samples
############################################################

## Physical params
params <- list(
  secpera = 31556926, # seconds per annum
  n = 3.0, # exponent in Glen's flow law
  rho_i = 917.0, # ice density
  rho_w = 1028.0, # sea water density
  g = 9.81 # gravity constant
  # A = 4.227e-25, #1.4579e-25, # flow rate parameter
)

params$m <- 1 / params$n
params$B <- 0.55 * 1e6 * params$secpera^params$m
params$A <- params$B^(-params$n)

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# 0. Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# 0. Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]

## Scaling units for friction coefficients
# fric_scale <- 1e6 * params$secpera^(1 / params$n)

## SMB data
smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)
params$as <- smb_avg # surface accumulation rate (m/s)

## Basal melt data
if (use_basal_melt_data) {
  melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
  # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
  avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
  melt_nonmissing <- which(!is.na(avg_melt_rate))
  avg_melt_rate[1:(melt_nonmissing[1] - 1)] <- -1 # impose a melt rate of 1 m/a upstream of the first non-missing value
  avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) # for the remaining part of the shelf just use the last non-missing value
  avg_melt_rate <- -avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
} else {
  avg_melt_rate <- rep(0, J)
}
params$ab <- avg_melt_rate # melt rate (m/s)

post_pred <- sim_obs(
  param_list = post_param_list,
  domain = domain,
  phys_params = params,
  years = 10, # sim_beds = T,
  warmup = 1,
  ini_thickness = ssa_steady$current_thickness,
  ini_velocity = ssa_steady$current_velocity,
  smb = smb_avg,
  basal_melt = avg_melt_rate
  # log_transform = log_transform
)

prior_pred <- sim_obs(
  param_list = prior_param_list,
  domain = domain,
  phys_params = params,
  years = 10, # sim_beds = T,
  warmup = 1,
  ini_thickness = ssa_steady$current_thickness,
  ini_velocity = ssa_steady$current_velocity,
  smb = smb_avg,
  basal_melt = avg_melt_rate
  # log_transform = log_transform
)

## Extract simulated observations
post_pred_out <- process_sim_results(sims = post_pred$results)
post_pred_obs <- post_pred_out$surface_obs_arr

prior_pred_out <- process_sim_results(sims = prior_pred$results)
prior_pred_obs <- prior_pred_out$surface_obs_arr

## Then average over the posterior samples
post_pred_obs_mean <- apply(post_pred_obs, c(2, 3, 4), mean)
prior_pred_obs_mean <- apply(prior_pred_obs, c(2, 3, 4), mean)

## Plot simulated surface observations
# s <- 2
png(filename = paste0(plot_dir, "post_pred_obs_", data_date, ".png"), width = 2000, height = 1500, res = 200)
par(mfrow = c(2, 2))
matplot(surf_elev_data,
  type = "l", lty = 1, col = "salmon",
  main = "Prior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
)
matlines(prior_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

matplot(velocity_data,
  type = "l", lty = 1, col = "salmon",
  main = "Prior predictive velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
)
matlines(prior_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

## Plot posterior predictive distribution with discrepancy adjustment
matplot(surf_elev_data,
  type = "l", lty = 1, col = "salmon",
  main = "Posterior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
)
matlines(post_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

matplot(velocity_data,
  type = "l", lty = 1, col = "salmon",
  main = "Posterior predictive surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
)
matlines(post_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)
dev.off()

## Model discrepancy
avg_vel_discr <- qread(file = paste0("./data/discrepancy/avg_vel_discr_", data_date, ".qs"))
avg_se_discr <- qread(file = paste0("./data/discrepancy/avg_se_discr_", data_date, ".qs"))

if (correct_model_discrepancy) {
  ## Add discrepancy to simulated observations
  post_pred_obs_adj <- post_pred_obs
  prior_pred_obs_adj <- prior_pred_obs
  for (s in 1:n_post_samples) {
    # s <- 1
    for (t in 1:dim(post_pred_obs)[3]) {
      post_pred_obs_adj[s, , t, 1] <- post_pred_obs_adj[s, , t, 1] + avg_se_discr
      post_pred_obs_adj[s, , t, 2] <- post_pred_obs_adj[s, , t, 2] + avg_vel_discr
      prior_pred_obs_adj[s, , t, 1] <- prior_pred_obs_adj[s, , t, 1] + avg_se_discr
      prior_pred_obs_adj[s, , t, 2] <- prior_pred_obs_adj[s, , t, 2] + avg_vel_discr
    }
  }

  ## Average over posterior samples
  adj_post_pred_obs_mean <- apply(post_pred_obs_adj, c(2, 3, 4), mean)
  adj_prior_pred_obs_mean <- apply(prior_pred_obs_adj, c(2, 3, 4), mean)

  ## Compare prior and posterior predictive distributions with discrepancy adjustment
  png(filename = paste0(plot_dir, "adj_post_pred_surface_obs_", data_date, ".png"), width = 2000, height = 1500, res = 200)

  ## Plot prior predictive distribution with discrepancy adjustment
  par(mfrow = c(2, 2))
  matplot(surf_elev_data,
    type = "l", lty = 1, col = "salmon",
    main = "Prior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
  )
  matlines(adj_prior_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
  legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

  matplot(velocity_data,
    type = "l", lty = 1, col = "salmon",
    main = "Prior predictive velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
  )
  matlines(adj_prior_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
  legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

  ## Plot posterior predictive distribution with discrepancy adjustment
  matplot(surf_elev_data,
    type = "l", lty = 1, col = "salmon",
    main = "Posterior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
  )
  matlines(adj_post_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
  legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

  matplot(velocity_data,
    type = "l", lty = 1, col = "salmon",
    main = "Posterior predictive surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
  )
  matlines(adj_post_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
  legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)
  dev.off()
}
