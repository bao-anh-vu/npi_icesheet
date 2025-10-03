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
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

data_date <- "20241111" # "20241103"
sets <- 1:50 # 6:20
# use_missing_pattern <- T
use_basal_melt_data <- T
correct_model_discrepancy <- F

setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

# if (use_missing_pattern) {
#   data_dir <- paste0("./data/training_data/", setsf, "/missing/")
#   output_dir <- paste0("./output/cnn/", setsf, "/missing/")
#   plot_dir <- paste0("./plots/cnn/", setsf, "/missing/")
# } else {
#   data_dir <- paste0("./data/training_data/", setsf, "/nonmissing/")
#   output_dir <- paste0("./output/cnn/", setsf, "/nonmissing/")
#   plot_dir <- paste0("./plots/cnn/", setsf, "/nonmissing/")
# }

data_dir <- paste0("./data/training_data/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")

## Directories
if (correct_model_discrepancy) {
    pred_output_dir <- paste0(output_dir, "pred/discr/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")
} else {
    pred_output_dir <- paste0(output_dir, "pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

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
params$B <- 0.6 * 1e6 * params$secpera^params$m
params$A <- params$B^(-params$n)

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# 0. Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# 0. Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
bed_obs_chosen <- bed_obs_df#[bed_obs_df$chosen == 1, ]

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

############################################################
##      Load posterior samples and compare with prior
############################################################

## Posterior samples
post_fric_samples <- qread(file = paste0(pred_output_dir, "fric_samples_real_", data_date, ".qs"))
post_bed_samples <- qread(file = paste0(pred_output_dir, "bed_samples_real_", data_date, ".qs"))

fric_post_mean <- qread(file = paste0(pred_output_dir, "pred_fric_real_", data_date, ".qs"))
bed_post_mean <- qread(file = paste0(pred_output_dir, "pred_bed_real_", data_date, ".qs"))

## Compare with simulations from the prior
# set <- sets[1]
# setf <- formatC(set, width = 2, flag = "0")
# prior_fric_samples <- t(qread(file = paste0("./data/training_data/friction_arr_", setf, "_", data_date, ".qs")))
# prior_bed_samples <- t(qread(file = paste0("./data/training_data/bed_arr_", setf, "_", data_date, ".qs")))
print("Simulating from prior...")
n_samples <- 100
bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
L <- t(chol(bed_prior$cov))
u_mat <- matrix(rnorm(nrow(L) * n_samples), nrow = nrow(L), ncol = n_samples) 
mean_mat <- matrix(rep(bed_prior$mean, n_samples), nrow = nrow(bed_prior$mean), ncol = n_samples)
prior_bed_samples <- mean_mat + L %*% u_mat 

prior_fric_samples <- simulate_friction2(
      nsim = n_samples, domain = domain) 

## Create parameter lists for posterior and prior samples

post_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = post_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = post_bed_samples[, r] #+ bed_mean
  )
})

prior_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = prior_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = prior_bed_samples[, r] #+ bed_mean
  )
})

## Plot prior vs posterior samples
png(file = paste0(plot_dir, "/posterior_vs_prior_", data_date, ".png"), width = 800, height = 800)
matplot(post_fric_samples[, 1:3], type = "l", col = "red", main = "Friction samples vs prior simulations", ylab = "Friction coeff (scaled)", xlab = "Sample index")
matlines(prior_fric_samples[, 1:3], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)

matplot(post_bed_samples[, 1:3], type = "l", col = "red", main = "Bed samples vs prior simulations", ylab = "Bedrock (m)", xlab = "Sample index")
matlines(prior_bed_samples[, 1:3], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)
dev.off()

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


warmup <- 0 # to let the model spin up
years <- dim(surf_elev_data)[2]  # simulate for the same number of years as the observations

post_pred <- sim_obs(
  param_list = post_param_list,
  domain = domain,
  phys_params = params,
  years = years, # sim_beds = T,
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
  years = years, # sim_beds = T,
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

# ## Combine simulations into posterior and prior predictive arrays
# post_pred_se <- post_pred_obs[, , , 1]
# post_pred_vel <- post_pred_obs[, , , 2]

# ## For each year, compute the mean and 95% CI of the surface elevation and velocity
# post_pred_se_mean <- apply(post_pred_se, c(2, 3), mean)

# post_pred_se_lci <- apply(post_pred_se, c(2, 3), quantile, probs = 0.025)
# post_pred_se_uci <- apply(post_pred_se, c(2, 3), quantile, probs = 0.975)

if (correct_model_discrepancy) {
  ## Model discrepancy
  avg_vel_discr <- qread(file = paste0("./data/discrepancy/avg_vel_discr_", data_date, ".qs"))
  avg_se_discr <- qread(file = paste0("./data/discrepancy/avg_se_discr_", data_date, ".qs"))

  ## Add discrepancy to simulated observations
  # post_pred_obs_adj <- post_pred_obs
  # prior_pred_obs_adj <- prior_pred_obs
  for (s in 1:n_samples) {
    # s <- 1
    for (t in 1:dim(post_pred_obs)[3]) {
      post_pred_obs[s, , t, 1] <- post_pred_obs[s, , t, 1] + avg_se_discr
      post_pred_obs[s, , t, 2] <- post_pred_obs[s, , t, 2] + avg_vel_discr
      prior_pred_obs[s, , t, 1] <- prior_pred_obs[s, , t, 1] + avg_se_discr
      prior_pred_obs[s, , t, 2] <- prior_pred_obs[s, , t, 2] + avg_vel_discr
    }
  }
}

########### Prior & Posterior predictive distributions ##########
prior_pred_obs_mean <- apply(prior_pred_obs, c(2, 3, 4), mean)
prior_pred_obs_lci <- apply(prior_pred_obs, c(2, 3, 4), quantile, probs = 0.025, na.rm = T)
prior_pred_obs_uci <- apply(prior_pred_obs, c(2, 3, 4), quantile, probs = 0.975, na.rm = T)

post_pred_obs_mean <- apply(post_pred_obs, c(2, 3, 4), mean)
post_pred_obs_lci <- apply(post_pred_obs, c(2, 3, 4), quantile, probs = 0.025, na.rm = T)
post_pred_obs_uci <- apply(post_pred_obs, c(2, 3, 4), quantile, probs = 0.975, na.rm = T)

## Turn into dataframe for plotting
### First for surface elevation

se_post_pred_df <- data.frame(
  domain = rep(domain, years),
  year = rep(1:years, each = J),
  obs = as.vector(surf_elev_data),
  prior_mean = as.vector(prior_pred_obs_mean[, , 1]),
  prior_lci = as.vector(prior_pred_obs_lci[, , 1]),
  prior_uci = as.vector(prior_pred_obs_uci[, , 1]),
  post_mean = as.vector(post_pred_obs_mean[, , 1]),
  post_lci = as.vector(post_pred_obs_lci[, , 1]),
  post_uci = as.vector(post_pred_obs_uci[, , 1])
)

## Then for surface velocity
vel_post_pred_df <- data.frame(
  domain = rep(domain, years),
  year = rep(1:years, each = J),
  obs = as.vector(velocity_data),
  prior_mean = as.vector(prior_pred_obs_mean[, , 2]),
  prior_lci = as.vector(prior_pred_obs_lci[, , 2]),
  prior_uci = as.vector(prior_pred_obs_uci[, , 2]),
  post_mean = as.vector(post_pred_obs_mean[, , 2]),
  post_lci = as.vector(post_pred_obs_lci[, , 2]),
  post_uci = as.vector(post_pred_obs_uci[, , 2])
)

## Now do facet plots of prior and posterior predictive distributions for surface elevation and velocity
library(ggplot2)
library(gridExtra)

png(filename = paste0(plot_dir, "post_pred_se_by_year.png"), width = 2000, height = 3000, res = 200)
p1 <- ggplot(se_post_pred_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.5) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  facet_wrap(~ year, ncol = 2) +
  xlim(0, 150e3) +
  labs(title = "Prior vs posterior predictive surface elevation", y = "Surface elevation (m)", x = "Domain") +
  theme_minimal()
print(p1)
dev.off()


## Same for the velocity
png(filename = paste0(plot_dir, "post_pred_vel_by_year.png"), width = 2000, height = 3000, res = 200)
p2 <- ggplot(vel_post_pred_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "red", alpha = 0.25) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  facet_wrap(~ year, ncol = 2, scales = "free_y") +
  xlim(0, 150e3) +
  # ylim(0, 4000) +
  labs(title = "Prior vs posterior predictive surface velocity", y = "Surface velocity (m/yr)", x = "Domain") +
  theme_minimal()
print(p2)
dev.off()

# browser()
# ## Plot simulated surface observations
# # s <- 1
# png(filename = paste0(plot_dir, "post_pred_obs_", data_date, ".png"), width = 2000, height = 1500, res = 200)
# par(mfrow = c(2, 2))
# matplot(surf_elev_data, ylim = c(0, 1500),
#   type = "l", lty = 1, col = "salmon",
#   main = "Prior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
# )
# matlines(prior_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
# # matlines(prior_pred_obs[s, , , 1], col = rgb(0, 0, 0, 0.3))
# legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

# matplot(velocity_data,
#   type = "l", lty = 1, col = "salmon",
#   main = "Prior predictive velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
# )
# matlines(prior_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
# # matlines(prior_pred_obs[s, , , 2], col = rgb(0, 0, 0, 0.3))

# legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

# ## Plot posterior predictive distribution with discrepancy adjustment
# matplot(surf_elev_data, ylim = c(0, 1500),
#   type = "l", lty = 1, col = "salmon",
#   main = "Posterior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
# )
# matlines(post_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
# # matlines(post_pred_obs[s, , , 1], col = rgb(0, 0, 0, 0.3))
# legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

# matplot(velocity_data,
#   type = "l", lty = 1, col = "salmon",
#   main = "Posterior predictive surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
# )
# matlines(post_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
# # matlines(post_pred_obs[s, , , 2], col = rgb(0, 0, 0, 0.3))
# legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)
# dev.off()

# if (correct_model_discrepancy) {
#   ## Model discrepancy
#   avg_vel_discr <- qread(file = paste0("./data/discrepancy/avg_vel_discr_", data_date, ".qs"))
#   avg_se_discr <- qread(file = paste0("./data/discrepancy/avg_se_discr_", data_date, ".qs"))

#   ## Add discrepancy to simulated observations
#   # post_pred_obs_adj <- post_pred_obs
#   # prior_pred_obs_adj <- prior_pred_obs
#   for (s in 1:n_samples) {
#     # s <- 1
#     for (t in 1:dim(post_pred_obs)[3]) {
#       post_pred_obs[s, , t, 1] <- post_pred_obs[s, , t, 1] + avg_se_discr
#       post_pred_obs[s, , t, 2] <- post_pred_obs[s, , t, 2] + avg_vel_discr
#       prior_pred_obs[s, , t, 1] <- prior_pred_obs[s, , t, 1] + avg_se_discr
#       prior_pred_obs[s, , t, 2] <- prior_pred_obs[s, , t, 2] + avg_vel_discr
#     }
#   }

#   ## Average over posterior samples
#   post_pred_obs_mean <- apply(post_pred_obs, c(2, 3, 4), mean)
#   prior_pred_obs_mean <- apply(prior_pred_obs, c(2, 3, 4), mean)

#   ## Compare prior and posterior predictive distributions with discrepancy adjustment
#   png(filename = paste0(plot_dir, "adj_post_pred_surface_obs_", data_date, ".png"), width = 2000, height = 1500, res = 200)

#   ## Plot prior predictive distribution with discrepancy adjustment
#   par(mfrow = c(2, 2))
#   matplot(surf_elev_data,
#     type = "l", lty = 1, col = "salmon",
#     main = "Prior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
#   )
#   matlines(prior_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
#   legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

#   matplot(velocity_data,
#     type = "l", lty = 1, col = "salmon",
#     main = "Prior predictive velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
#   )
#   matlines(prior_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
#   legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

#   ## Plot posterior predictive distribution with discrepancy adjustment
#   matplot(surf_elev_data,
#     type = "l", lty = 1, col = "salmon",
#     main = "Posterior predictive surface elevation", ylab = "Surface elevation (m)", xlab = "Domain"
#   )
#   matlines(post_pred_obs_mean[, , 1], col = rgb(0, 0, 0, 0.3))
#   legend("topright", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)

#   matplot(velocity_data,
#     type = "l", lty = 1, col = "salmon",
#     main = "Posterior predictive surface velocity", ylab = "Surface velocity (m/yr)", xlab = "Domain"
#   )
#   matlines(post_pred_obs_mean[, , 2], col = rgb(0, 0, 0, 0.3))
#   legend("topleft", legend = c("Simulated", "Observed"), col = c(rgb(0, 0, 0, 0.3), "salmon"), lty = 1)
#   dev.off()
# }

## Compute RMSE
rmse <- function(obs, sim) {
  sqrt(mean((obs - sim)^2, na.rm = T))
}

prior_se_rmse <- rmse(se_post_pred_df$obs, se_post_pred_df$prior_mean)
prior_vel_rmse <- rmse(vel_post_pred_df$obs, vel_post_pred_df$prior_mean)
post_se_rmse <- rmse(se_post_pred_df$obs, se_post_pred_df$post_mean)
post_vel_rmse <- rmse(vel_post_pred_df$obs, vel_post_pred_df$post_mean)

## Collect RMSE into a data frame
rmse_df <- data.frame(
  Pred = c("Prior", "Posterior"),
  surface_elev_RMSE = c(prior_se_rmse, post_se_rmse),
  velocity_RMSE = c(prior_vel_rmse, post_vel_rmse)
)

## Save as a .csv file
write.csv(rmse_df, file = paste0(plot_dir, "rmse_summary_", data_date, ".csv"), row.names = F)
