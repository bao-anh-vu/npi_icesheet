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

data_date <- "20241111" #"20241103"
sets <- 51:60 #6:20
use_missing_pattern <- T

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
fric_samples <- qread(file = paste0(output_dir, "/fric_samples_real_", data_date, ".qs"))
bed_samples <- qread(file = paste0(output_dir, "/bed_samples_real_", data_date, ".qs"))

fric_post_mean <- qread(file = paste0(output_dir, "/pred_fric_real_", data_date, ".qs"))
bed_post_mean <- qread(file = paste0(output_dir, "/pred_bed_real_", data_date, ".qs"))
    

## Compare with simulations from the prior
set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")
fric_sims <- qread(file = paste0("./data/training_data/friction_arr_", setf, "_", data_date, ".qs"))


png(file = paste0(plot_dir, "/posterior_vs_prior_", data_date, ".png"), width = 800, height = 800)
matplot(fric_samples[,1:3], type = 'l', col = 'red', main = "Friction samples vs prior simulations", ylab = "Friction coeff (scaled)", xlab = "Sample index")
matlines(t(fric_sims)[, 1:3], col = 'blue')
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)
# plot(bed_samples[,1], type = 'l', col = 'blue', main = "Bed samples vs prior simulations", ylab = "Bed coeff (scaled)", xlab = "Sample index")
# lines(bed_sims[1, ], col = 'red')
dev.off()

param_list <- lapply(1:3, function(r) {
    list(
    friction = fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = bed_samples[, r] #+ bed_mean
    )
})

param_mean <- list()
param_mean[[1]] <- list(
    friction = as.vector(fric_post_mean), #* 1e6 * params$secpera^(1 / params$n),
    bedrock = as.vector(bed_post_mean) #+ bed_mean
)

## Plot parameters
png(filename = paste0(plot_dir, "/fric_bed_samples_real_", data_date, ".png"), width = 800, height = 400)
par(mfrow = c(2, 1))
plot(param_list[[1]]$friction, type = "l", main = "Friction sample", ylab = "Friction coeff (scaled)", xlab = "Sample index")
plot(param_list[[1]]$bedrock, type = "l", main = "Bed sample", ylab = "Friction coeff (scaled)", xlab = "Sample index")
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
params$B <- 1.4 * 1e6 * params$secpera^params$m
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

avg_melt_rate <- rep(0, length(domain)) # assume no melt for now
params$ab <- avg_melt_rate # basal melt rate (m/s)

sim_results <- sim_obs(
        param_list = param_mean, #param_list,
        domain = domain,
        phys_params = params,
        years = 10, # sim_beds = T,
        # warmup = warmup,
        ini_thickness = ssa_steady$current_thickness,
        ini_velocity = ssa_steady$current_velocity,
        smb = smb_avg,
        basal_melt = avg_melt_rate
        # log_transform = log_transform
      )

## Extract simulated observations
generated_data <- process_sim_results(sims = sim_results$results)
surface_obs_arr <- generated_data$surface_obs_arr

## Plot simulated surface observations
s <- 1
png(filename = paste0(plot_dir, "/predictive_surface_obs_", data_date, ".png"), width = 800, height = 800)
par(mfrow = c(2, 1))
matplot(surface_obs_arr[s,,,1], type = 'l', lty = 1, col = rgb(0,0,0,0.3), 
main = "Simulated surface elevation", ylab = "Surface elevation (m)", xlab = "Domain")
matlines(surf_elev_data, col = 'salmon')
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0,0,0,0.3), 'salmon'), lty = 1)

matplot(surface_obs_arr[s,,,2], type = 'l', lty = 1, col = rgb(0,0,0,0.3),
main = "Simulated surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain")
matlines(velocity_data, col = 'salmon')
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0,0,0,0.3), 'salmon'), lty = 1)
dev.off()

## Model discrepancy
avg_vel_discr <- qread(file = paste0("./data/discrepancy/avg_vel_discr_", data_date, ".qs"))
avg_se_discr <- qread(file = paste0("./data/discrepancy/avg_se_discr_", data_date, ".qs"))

## Add discrepancy to simulated observations
surface_obs_arr_adj <- surface_obs_arr
# for (s in 1:length(param_list)) {
s <- 1
    for (t in 1:dim(surface_obs_arr)[3]) {
        surface_obs_arr_adj[s, , t, 1] <- surface_obs_arr[s, , t, 1] + avg_se_discr
        surface_obs_arr_adj[s, , t, 2] <- surface_obs_arr[s, , t, 2] + avg_vel_discr
    }
# }

png(filename = paste0(plot_dir, "/predictive_adj_surface_obs_", data_date, ".png"), width = 800, height = 800)
par(mfrow = c(2, 1))
matplot(surface_obs_arr_adj[s,,,1], type = 'l', lty = 1, col = rgb(0,0,0,0.3), 
main = "Simulated surface elevation", ylab = "Surface elevation (m)", xlab = "Domain")
matlines(surf_elev_data, col = 'salmon')
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0,0,0,0.3), 'salmon'), lty = 1)

matplot(surface_obs_arr[s,,,2], type = 'l', lty = 1, col = rgb(0,0,0,0.3),
main = "Simulated surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain")
matlines(velocity_data, col = 'salmon')
legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0,0,0,0.3), 'salmon'), lty = 1)
dev.off()
