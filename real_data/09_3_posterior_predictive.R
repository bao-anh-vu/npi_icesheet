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
library(FRK)
library(qs)
library(parallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

source("./source/sim_obs.R")
# source("./source/cond_sim_gp.R")
source("./source/get_ini_thickness.R")
source("./source/process_sim_results.R")
source("./source/fit_basis.R")
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
source("./source/comparison_metrics.R")

data_date <- "20241111" # "20241103"
sets <- 101:150 # 6:20
# use_missing_pattern <- T
use_basal_melt_data <- T
correct_model_discrepancy <- T

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
# velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

## Physical params
# params <- list(
#   secpera = 31556926, # seconds per annum
#   n = 3.0, # exponent in Glen's flow law
#   rho_i = 917.0, # ice density
#   rho_w = 1028.0, # sea water density
#   g = 9.81 # gravity constant
#   # A = 4.227e-25, #1.4579e-25, # flow rate parameter
# )

# params$m <- 1 / params$n
# params$B <- 0.6 * 1e6 * params$secpera^params$m
# params$A <- params$B^(-params$n)

# ## SMB data
# smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
# smb_avg <- colMeans(smb_data_racmo, na.rm = T)
# params$as <- smb_avg # surface accumulation rate (m/s)

# ## Basal melt data
# if (use_basal_melt_data) {
#   melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
#   # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
#   avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
#   melt_nonmissing <- which(!is.na(avg_melt_rate))
#   avg_melt_rate[1:(melt_nonmissing[1] - 1)] <- -1 # impose a melt rate of 1 m/a upstream of the first non-missing value
#   avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) # for the remaining part of the shelf just use the last non-missing value
#   avg_melt_rate <- -avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
# } else {
#   avg_melt_rate <- rep(0, J)
# }
# params$ab <- avg_melt_rate # melt rate (m/s)

params <- qread(file = paste0("./data/training_data/", "/phys_params_", data_date, ".qs"))

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# 0. Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# 0. Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
# bed_obs_chosen <- bed_obs_df#[bed_obs_df$chosen == 1, ]

## BedMachine data for comparison
bedmachine <- qread("./data/bedmachine/flowline_bedmachine.qs")

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
set.seed(2025)
n_samples <- 1000
bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
L <- t(chol(bed_prior$cov))
u_mat <- matrix(rnorm(nrow(L) * n_samples), nrow = nrow(L), ncol = n_samples) 
mean_mat <- matrix(rep(bed_prior$mean, n_samples), nrow = nrow(bed_prior$mean), ncol = n_samples)
prior_bed_samples <- mean_mat + L %*% u_mat 

prior_fric_samples <- simulate_friction2(
      nsim = n_samples, domain = domain) 

## Fit basis functions to prior samples
n_fric_basis <- 120
n_bed_basis <- 150
prior_fric_basis <- fit_friction_basis(
        nbasis = n_fric_basis,
        domain = domain,
        fric_arr = t(prior_fric_samples),
        log_transform = T,
        lengthscale = 3e3
    )

bed_arr_demean <- prior_bed_samples - mean_mat

## Fit basis to de-meaned bedrock
prior_bed_basis <- fit_bed_basis(
    nbasis = n_bed_basis, domain = domain,
    bed_arr = t(bed_arr_demean),
    lengthscale = 2.5e3
)
prior_bed_basis$mean <- bed_prior$mean

prior_fric_samples_basis <- t(exp(prior_fric_basis$fitted_values))
prior_bed_samples_basis <- t(prior_bed_basis$fitted_values) + mean_mat

## Create parameter lists for posterior and prior samples

post_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = post_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = post_bed_samples[, r]#, #+ bed_mean
    # ini_thickness = ssa_steady$current_thickness,
    # ini_velocity = ssa_steady$current_velocity
  )
})

prior_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = prior_fric_samples_basis[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = prior_bed_samples_basis[, r]#,
    # ini_thickness = ssa_steady$current_thickness,
    # ini_velocity = ssa_steady$current_velocity
  )
})

## Plot prior vs posterior samples
png(file = paste0(plot_dir, "posterior_vs_prior_", data_date, ".png"), width = 800, height = 800)
par(mfrow = c(2, 1))
sims <- 1:3
matplot(post_fric_samples[, sims], type = "l", col = "red", main = "Friction samples vs prior simulations", ylab = "Friction coeff (scaled)", xlab = "Sample index")
matlines(prior_fric_samples_basis[, sims], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)

matplot(post_bed_samples[, sims], type = "l", col = "red", main = "Bed samples vs prior simulations", ylab = "Bedrock (m)", xlab = "Sample index")
matlines(prior_bed_samples_basis[, sims], col = "blue")
legend("topright", legend = c("Posterior samples", "Prior simulations"), col = c("red", "blue"), lty = 1)
dev.off()

## Same plot but in ggplot

n_plot_samples <- 10
prior_fric_df <- data.frame(x = rep(domain/1000, n_plot_samples),
                      sample = rep(1:n_plot_samples, each = length(domain)),
                      fric = as.vector(prior_fric_samples_basis[, 1:n_plot_samples]))
post_fric_df <- data.frame(x = rep(domain/1000, n_plot_samples),
                      sample = rep(1:n_plot_samples, each = length(domain)),
                      fric = as.vector(post_fric_samples[, 1:n_plot_samples]))
# post_fric_df_sample <- post_fric_df %>% filter(sample %in% 1:n_plot_samples)
# post_fric_df_wide <- pivot_wider(post_fric_df, names_from = sample, values_from = fric, names_prefix = "sample")
# post_fric_df_long <- pivot_longer(post_fric_df_wide, cols = starts_with("sample"), names_to = "sample", values_to = "fric")
post_fric_plot <- ggplot() +
  geom_line(data = prior_fric_df, aes(x = x, y = fric, group = sample), col = "lightblue") +
  geom_line(data = post_fric_df, aes(x = x, y = fric, group = sample), col = "salmon") +
  geom_line(alpha = 0.3) +
  theme_bw() +
  xlim(c(0, 150)) +
  xlab("Dist. along flowline (km)") +
  ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")")) +
  theme(text = element_text(size = 24))

# png(filename = paste0(plot_dir, "fric_samples_real_gg_", data_date, ".png"), width = 1000, height = 600, res = 150)
# print(post_fric_plot)
# dev.off()

prior_bed_df <- data.frame(x = rep(domain/1000, n_plot_samples),
                      sample = rep(1:n_plot_samples, each = length(domain)),
                      bed = as.vector(prior_bed_samples_basis[, 1:n_plot_samples]))

post_bed_df <- data.frame(x = rep(domain/1000, n_plot_samples),
                      sample = rep(1:n_plot_samples, each = length(domain)),
                      bed = as.vector(post_bed_samples[, 1:n_plot_samples]))

bedmachine_df <- data.frame(x = domain/1000, bed = bedmachine$bed_avg)

# post_bed_df_sample <- post_bed_df %>% filter(sample %in% 1:n_plot_samples)

post_bed_plot <- ggplot() +
  geom_line(data = prior_bed_df, aes(x = x, y = bed, group = sample), col = "lightblue") +
  geom_line(data = post_bed_df, aes(x = x, y = bed, group = sample), col = "salmon") +
  geom_line(data = bedmachine_df, aes(x = x, y = bed), col = "blue", lwd = 1, lty = 2) +
  geom_point(data = bed_obs_df, aes(x = loc/1000, y = bed_elev), col = "black", size = 2) +
  theme_bw() +
  xlim(c(0, 150)) +
  ylim(c(-1500, -500)) +
  xlab("Dist. along flowline (km)") +
  ylab("Bed (m)") +
  theme(text = element_text(size = 24))

png(filename = paste0(plot_dir, "param_samples_real_gg_", data_date, ".png"), width = 1000, height = 1000, res = 150)
grid.arrange(grobs = list(post_bed_plot, post_fric_plot), nrow = 2, ncol = 1)
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
  warmup = warmup,
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
  warmup = warmup,
  ini_thickness = ssa_steady$current_thickness,
  ini_velocity = ssa_steady$current_velocity,
  smb = smb_avg,
  basal_melt = avg_melt_rate
  # log_transform = log_transform
)

prior_bad_sims <- prior_pred$bad_sims
if (length(prior_bad_sims) > 0) {
    cat(length(prior_bad_sims), " prior simulations failed \n")
    prior_pred$results <- prior_pred$results[-prior_bad_sims]
}

post_bad_sims <- post_pred$bad_sims
if (length(post_bad_sims) > 0) {
    cat(length(post_bad_sims), " posterior simulations failed \n")
    post_pred$results <- post_pred$results[-post_bad_sims]
}

## Extract simulated observations
prior_pred_out <- process_sim_results(sims = prior_pred$results)
prior_pred_obs <- prior_pred_out$surface_obs_arr

post_pred_out <- process_sim_results(sims = post_pred$results)
post_pred_obs <- post_pred_out$surface_obs_arr


# ## Combine simulations into posterior and prior predictive arrays
# post_pred_se <- post_pred_obs[, , , 1]
# post_pred_vel <- post_pred_obs[, , , 2]

# ## For each year, compute the mean and 95% CI of the surface elevation and velocity
# post_pred_se_mean <- apply(post_pred_se, c(2, 3), mean)

# post_pred_se_lci <- apply(post_pred_se, c(2, 3), quantile, probs = 0.025)
# post_pred_se_uci <- apply(post_pred_se, c(2, 3), quantile, probs = 0.975)

if (correct_model_discrepancy) {

  if (avg_over_time) {
    file_tag <- "avg_"
  } else {
      file_tag <- ""
  }

  ## Model discrepancy
  vel_discr <- qread(file = paste0("./data/discrepancy/", file_tag, "vel_discr_", data_date, ".qs"))
  se_discr <- qread(file = paste0("./data/discrepancy/", file_tag, "se_discr_", data_date, ".qs"))

  # qsave(adj_se_mat, file = paste0(data_dir, "surface_elev/adj_se_mat_", file_tag, data_date, ".qs"))
  # qsave(adj_vel_mat, file = paste0(data_dir, "velocity/adj_vel_mat_", file_tag, data_date, ".qs"))

  ## Add discrepancy to simulated observations
  for (s in 1:length(prior_pred$results)) {
      prior_pred_obs[s, , , 1] <- prior_pred_obs[s, , t, 1] + se_discr
      prior_pred_obs[s, , , 2] <- prior_pred_obs[s, , t, 2] + vel_discr
  }
  
  for (s in 1:length(post_pred$results)) {
    # s <- 1
      post_pred_obs[s, , , 1] <- post_pred_obs[s, , , 1] + se_discr
      post_pred_obs[s, , , 2] <- post_pred_obs[s, , , 2] + vel_discr
  }

  # for (s in 1:length(prior_pred$results)) {
  #   # s <- 1
  #   for (t in 1:dim(prior_pred_obs)[3]) {
  #     prior_pred_obs[s, , t, 1] <- prior_pred_obs[s, , t, 1] + se_discr
  #     prior_pred_obs[s, , t, 2] <- prior_pred_obs[s, , t, 2] + vel_discr
  #   }
  # }
  
  # for (s in 1:length(post_pred$results)) {
  #   # s <- 1
  #   for (t in 1:dim(post_pred_obs)[3]) {
  #     post_pred_obs[s, , t, 1] <- post_pred_obs[s, , t, 1] + se_discr
  #     post_pred_obs[s, , t, 2] <- post_pred_obs[s, , t, 2] + vel_discr
  #   }
  # }

}

## Save simulated observations
qsave(prior_pred_obs, file = paste0(pred_output_dir, "prior_pred_obs_", data_date, ".qs"))
qsave(post_pred_obs, file = paste0(pred_output_dir, "post_pred_obs_", data_date, ".qs"))

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
  domain = rep(domain/1000, years),
  year = rep(2010 + 0:(years-1), each = J),
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
  domain = rep(domain/1000, years),
  year = rep(2010 + 0:(years-1), each = J),
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
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.4) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  facet_wrap(~ year, ncol = 2) +
  xlim(0, 150) +
  labs(title = "Prior vs posterior predictive surface elevation by year", 
      y = "Surface elevation (m)", x = "Dist. along flowline (km)") +
  theme_bw() + 
  theme(text = element_text(size = 20))
print(p1)
dev.off()

## Do the same plot, but save each year individually
yearnum <- 2010:2020
for (yr in yearnum) {

  se_df <- se_post_pred_df[se_post_pred_df$year == yr, ]

  png(filename = paste0(plot_dir, "post_pred_se_year", formatC(yr, width = 2, flag = "0"), "_", data_date, ".png"), width = 1000, height = 600, res = 150)

  p_yr <- ggplot(se_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.5) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  # facet_wrap(~ year, ncol = 2) +
  xlim(0, 150) +
  labs(title = paste0("Year ", yr), y = "Surface elevation (m)", x = "Dist. along flowline (km)") +
  theme_bw() +
  theme(text = element_text(size = 20))
  
  print(p_yr)

  dev.off()
}


## Same for the velocity
png(filename = paste0(plot_dir, "post_pred_vel_by_year.png"), width = 2000, height = 3000, res = 200)
p2 <- ggplot(vel_post_pred_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.4) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  facet_wrap(~ year, ncol = 2, scales = "free_y") +
  xlim(0, 150) +
  # ylim(0, 4000) +
  labs(title = "Prior vs posterior predictive surface velocity", 
      y = "Surface velocity (m/yr)", x = "Dist. along flowline (km)") +
  theme_bw() +
  theme(text = element_text(size = 20))
print(p2)
dev.off()

## Individual plots for each year
for (yr in yearnum) {

  vel_df <- vel_post_pred_df[vel_post_pred_df$year == yr, ]

  png(filename = paste0(plot_dir, "post_pred_vel_year", formatC(yr, width = 2, flag = "0"), "_", data_date, ".png"), width = 1000, height = 600, res = 150)

  p_yr <- ggplot(vel_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.5) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  # facet_wrap(~ year, ncol = 2) +
  xlim(0, 150) +
  labs(title = paste0("Year ", yr), y = "Velocity (m/yr)", x = "Dist. along flowline (km)") +
  theme_bw() +
  theme(text = element_text(size = 20))
  
  print(p_yr)

  dev.off()
}

## Compute RMSE
# rmse <- function(obs, sim) {
#   sqrt(mean((obs - sim)^2, na.rm = T))
# }

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


## Compute CRPS

### Extract all final-year surface elevation obs from posterior predictive simulations
prior_pred_se_final <- prior_pred_obs[, , years, 1]
post_pred_se_final <- post_pred_obs[, , years, 1] 
obs_se_final <- surf_elev_data[, years]

### Discard observation locations with missing data
valid_locs <- which(!is.na(obs_se_final))
prior_pred_se_final <- prior_pred_se_final[, valid_locs]
post_pred_se_final <- post_pred_se_final[, valid_locs]
obs_se_final <- obs_se_final[valid_locs]

prior_crps_se_fin <- crps_multivariate(prior_pred_se_final,
                                obs_se_final)
post_crps_se_fin <- crps_multivariate(post_pred_se_final,
                              obs_se_final)

## Now calculate CRPS for velocity
prior_pred_vel_final <- prior_pred_obs[, , years, 2]
post_pred_vel_final <- post_pred_obs[, , years, 2]
obs_vel_final <- velocity_data[, years]

### Discard observation locations with missing data
valid_locs_vel <- valid_locs #which(!is.na(obs_vel_final))
prior_pred_vel_final <- prior_pred_vel_final[, valid_locs_vel]
post_pred_vel_final <- post_pred_vel_final[, valid_locs_vel]
obs_vel_final <- obs_vel_final[valid_locs_vel]

prior_crps_vel_fin <- crps_multivariate(prior_pred_vel_final,
                                obs_vel_final)
post_crps_vel_fin <- crps_multivariate(post_pred_vel_final,
                              obs_vel_final)

## Collect CRPS into a data frame
crps_df <- data.frame(
  Pred = c("Prior", "Posterior"),
  surface_elev_CRPS = c(prior_crps_se_fin, post_crps_se_fin),
  velocity_CRPS = c(prior_crps_vel_fin, post_crps_vel_fin)
)

## Save as a .csv file
write.csv(crps_df, file = paste0(plot_dir, "crps_summary_", data_date, ".csv"), row.names = F)
