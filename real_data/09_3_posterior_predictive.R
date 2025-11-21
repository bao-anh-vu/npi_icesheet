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
library(dplyr)
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

## Flags/settings
data_date <- "20241111" # "20241103"
sets <- 51:100 # 6:20
# use_missing_pattern <- T
use_basal_melt_data <- T
correct_model_discrepancy <- T
correct_velocity_discrepancy <- T # can correct velocity discrepancy as well, if F will just correct surface elevation discrepancy
avg_over_time <- T
warmup <- 5
n_samples <- 100

## Directories
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./data/training_data/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")

if (correct_model_discrepancy) {
  if (avg_over_time) {
    pred_output_dir <- paste0(output_dir, "pred/discr_avg/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr_avg/")
  } else {
    pred_output_dir <- paste0(output_dir, "pred/discr/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")
  }
} else {
  pred_output_dir <- paste0(output_dir, "pred/")
  plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
# velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

## Measurement error standard deviations
vel_err_sd <- qread(file = paste0("./data/velocity/vel_err_sd_", data_date, ".qs"))

# Load physical parameters (SMB, basal melt rate, etc.)
params <- qread(file = paste0("./data/training_data/", "/phys_params_", data_date, ".qs"))

# Load ice sheet at steady state
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Observed GL position
gl_pos <- qread(file = paste0("data/grounding_line/gl_pos.qs"))
gl_obs <- domain[gl_pos$ind]

# Read bed observations
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

# fric_post_mean <- qread(file = paste0(pred_output_dir, "pred_fric_real_", data_date, ".qs"))
# bed_post_mean <- qread(file = paste0(pred_output_dir, "pred_bed_real_", data_date, ".qs"))

## Compare with simulations from the prior
# set <- sets[1]
# setf <- formatC(set, width = 2, flag = "0")
# prior_fric_samples <- t(qread(file = paste0("./data/training_data/friction_arr_", setf, "_", data_date, ".qs")))
# prior_bed_samples <- t(qread(file = paste0("./data/training_data/bed_arr_", setf, "_", data_date, ".qs")))
print("Simulating from prior...")
set.seed(2025)

bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
L <- t(chol(bed_prior$cov))
u_mat <- matrix(rnorm(nrow(L) * n_samples), nrow = nrow(L), ncol = n_samples)
mean_mat <- matrix(rep(bed_prior$mean, n_samples), nrow = nrow(bed_prior$mean), ncol = n_samples)
prior_bed_samples <- mean_mat + L %*% u_mat

prior_fric_samples <- simulate_friction2(
  nsim = n_samples, domain = domain
)

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

qsave(prior_fric_samples_basis, file = paste0(pred_output_dir, "prior_fric_samples_basis_", data_date, ".qs"))
qsave(prior_bed_samples_basis, file = paste0(pred_output_dir, "prior_bed_samples_basis_", data_date, ".qs"))

# prior_fric_samples_basis <- qread(file = paste0(pred_output_dir, "prior_fric_samples_basis_", data_date, ".qs"))
# prior_bed_samples_basis <- qread(file = paste0(pred_output_dir, "prior_bed_samples_basis_", data_date, ".qs"))

## Create parameter lists for posterior and prior samples

post_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = post_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = post_bed_samples[, r] # , #+ bed_mean
    # ini_thickness = ssa_steady$current_thickness,
    # ini_velocity = ssa_steady$current_velocity
  )
})

prior_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = prior_fric_samples_basis[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = prior_bed_samples_basis[, r] # ,
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

n_plot_samples <- 4
prior_fric_df <- data.frame(
  x = rep(domain / 1000, n_plot_samples),
  sample = rep(1:n_plot_samples, each = length(domain)),
  fric = as.vector(prior_fric_samples_basis[, 1:n_plot_samples])
)
post_fric_df <- data.frame(
  x = rep(domain / 1000, n_plot_samples),
  sample = rep(1:n_plot_samples, each = length(domain)),
  fric = as.vector(post_fric_samples[, 1:n_plot_samples])
)
# post_fric_df_sample <- post_fric_df %>% filter(sample %in% 1:n_plot_samples)
# post_fric_df_wide <- pivot_wider(post_fric_df, names_from = sample, values_from = fric, names_prefix = "sample")
# post_fric_df_long <- pivot_longer(post_fric_df_wide, cols = starts_with("sample"), names_to = "sample", values_to = "fric")
post_fric_plot <- ggplot() +
  geom_line(data = prior_fric_df, aes(x = x, y = fric, group = sample), col = "lightblue") +
  geom_line(data = post_fric_df, aes(x = x, y = fric, group = sample), col = "salmon") +
  geom_line(alpha = 0.3) +
  theme_bw() +
  xlim(c(0, 150)) +
  xlab("Distance along flowline (km)") +
  ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")")) +
  theme(text = element_text(size = 24))

# png(filename = paste0(plot_dir, "fric_samples_real_gg_", data_date, ".png"), width = 1000, height = 600, res = 150)
# print(post_fric_plot)
# dev.off()

prior_bed_df <- data.frame(
  x = rep(domain / 1000, n_plot_samples),
  sample = rep(1:n_plot_samples, each = length(domain)),
  bed = as.vector(prior_bed_samples_basis[, 1:n_plot_samples])
)

post_bed_df <- data.frame(
  x = rep(domain / 1000, n_plot_samples),
  sample = rep(1:n_plot_samples, each = length(domain)),
  bed = as.vector(post_bed_samples[, 1:n_plot_samples])
)

bedmachine_df <- data.frame(x = domain / 1000, bed = bedmachine$bed_avg)

# post_bed_df_sample <- post_bed_df %>% filter(sample %in% 1:n_plot_samples)

post_bed_plot <- ggplot() +
  geom_line(data = prior_bed_df, aes(x = x, y = bed, group = sample), col = "lightblue") +
  geom_line(data = post_bed_df, aes(x = x, y = bed, group = sample), col = "salmon") +
  geom_line(data = bedmachine_df, aes(x = x, y = bed), col = "blue", lwd = 1, lty = 2) +
  geom_point(data = bed_obs_df, aes(x = loc / 1000, y = bed_elev), col = "black", size = 2) +
  theme_bw() +
  xlim(c(0, 150)) +
  ylim(c(-1500, -500)) +
  xlab("Distance along flowline (km)") +
  ylab("Bed (m)") +
  theme(text = element_text(size = 24))

png(filename = paste0(plot_dir, "param_samples_real_gg_", data_date, ".png"), width = 1000, height = 1000, res = 150)
grid.arrange(grobs = list(post_bed_plot, post_fric_plot), nrow = 2, ncol = 1)
dev.off()

############################################################
##      Simulate observations based on posterior samples
############################################################

years <- dim(surf_elev_data)[2] # simulate for the same number of years as the observations

post_pred <- sim_obs(
  param_list = post_param_list,
  domain = domain,
  phys_params = params,
  years = years, # sim_beds = T,
  warmup = warmup,
  # ini_thickness = ssa_steady$current_thickness,
  ini_surface = ssa_steady$current_top_surface,
  ini_velocity = ssa_steady$current_velocity,
  vel_err_sd = vel_err_sd
  # smb = smb_avg,
  # basal_melt = avg_melt_rate
  # log_transform = log_transform
)

prior_pred <- sim_obs(
  param_list = prior_param_list,
  domain = domain,
  phys_params = params,
  years = years, # sim_beds = T,
  warmup = warmup,
  # ini_thickness = ssa_steady$current_thickness,
  ini_surface = ssa_steady$current_top_surface,
  ini_velocity = ssa_steady$current_velocity,
  vel_err_sd = vel_err_sd
  # smb = smb_avg,
  # basal_melt = avg_melt_rate
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

qsave(prior_pred_obs, file = paste0(pred_output_dir, "prior_pred_obs_unadj_", data_date, ".qs"))

## Model discrepancy
if (correct_model_discrepancy) {
  # avg_vel_discr <- rowMeans(vel_discr, na.rm = T)
  # avg_se_discr <- rowMeans(se_discr, na.rm = T)

  if (avg_over_time) { # just use average discrepancy for all years
    vel_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/vel_discr_avg_", data_date, ".qs"))
    se_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/se_discr_avg_", data_date, ".qs"))

    # ## Fit a spline through the avg discrepancy
    #   library(mgcv)
    #   spline_fit <- gam(avg_se_discr ~ s(domain, k = 20))
    #   se_discr_smooth <- predict(spline_fit, newdata = data.frame(domain = domain))
    #   # se_discr_smooth <- se_discr_smooth[!is.na(avg_se_discr)]

    #   spline_fit <- gam(avg_vel_discr ~ s(domain, k = 20))
    #   vel_discr_smooth <- predict(spline_fit, newdata = data.frame(domain = domain))

    #   se_discr_smooth[is.na(avg_se_discr)] <- 0
    #   vel_discr_smooth[is.na(avg_vel_discr)] <- 0

    #   vel_discr_mat <- matrix(rep(vel_discr_smooth, years), nrow = length(vel_discr_smooth), ncol = years)
    #   se_discr_mat <- matrix(rep(se_discr_smooth, years), nrow = length(se_discr_smooth), ncol = years)

    # vel_discr_mat <- matrix(rep(avg_vel_discr, years), nrow = length(avg_vel_discr), ncol = years)
    # se_discr_mat <- matrix(rep(avg_se_discr, years), nrow = length(avg_se_discr), ncol = years)
  } else {
    # use yearly discrepancy, except for the one year that was left out for validation
    # for that year, fit a regression model for each spatial location based on the previous years' discrepancy

    # se_discr_mat <- cbind(se_discr, se_discr[, years-1])

    ## Do facet plot of velocity discrepancy for selected spatial locations
    locs <- seq(from = 1, to = 1000, by = 100)
    vel_discr_df <- data.frame(
      year = rep(1:(years - 1), each = length(locs)),
      location = rep(locs, times = years - 1),
      vel_discr = as.vector(vel_discr[locs, ])
    )

    vel_discr_plot <- ggplot(data = vel_discr_df, aes(x = year, y = vel_discr)) +
      geom_line() +
      facet_wrap(~location, ncol = 2) +
      theme_bw() +
      labs(x = "Year", y = "Velocity discrepancy (m/yr)") +
      ggtitle("Velocity discrepancy at selected grid pts") +
      theme(text = element_text(size = 18))
    png(paste0(plot_dir, "vel_discrepancy_facets_", data_date, ".png"), width = 1000, height = 1500, res = 150)
    print(vel_discr_plot)
    dev.off()

    ## Same plot for the surface elevation discrepancy
    locs_se <- seq(from = 1, to = 1000, by = 100)
    se_discr_df <- data.frame(
      year = rep(1:(years - 1), each = length(locs_se)),
      location = rep(locs_se, times = years - 1),
      se_discr = as.vector(se_discr[locs_se, ])
    )
    se_discr_plot <- ggplot(data = se_discr_df, aes(x = year, y = se_discr)) +
      geom_line() +
      facet_wrap(~location, ncol = 2) +
      theme_bw() +
      labs(x = "Year", y = "Surface elevation discrepancy (m)") +
      ggtitle("Surface elevation discrepancy at selected grid pts") +
      theme(text = element_text(size = 18))
    png(paste0(plot_dir, "se_discrepancy_facets_", data_date, ".png"), width = 1000, height = 1500, res = 150)
    print(se_discr_plot)
    dev.off()

    se_discr_fin <- rep(NA, nrow(se_discr))
    non_na_locs <- which(!is.na(se_discr[, 1]))
    for (j in non_na_locs) {
      # j <- 1
      # vel_lm <- lm(vel_discr[j, 1:(years-1)] ~ I(2010 + 0:(years-2)))
      # vel_pred <- predict(vel_lm, newdata = data.frame(year = 2010 + (years - 1)))
      # vel_discr_mat[j, years] <- vel_pred

      se_discr_df <- data.frame(year = 1:ncol(se_discr), se = se_discr[j, ])
      se_lm <- lm(se ~ year, data = na.omit(se_discr_df))
      se_pred <- predict(se_lm, newdata = data.frame(year = years))
      se_discr_fin[j] <- se_pred
    }
    se_discr_mat <- cbind(se_discr, se_discr_fin)

    ## do the same for the velocity
    vel_discr_fin <- rep(NA, nrow(vel_discr))
    non_na_locs_vel <- which(!is.na(vel_discr[, ncol(vel_discr)]))

    for (j in non_na_locs_vel) {
      vel_discr_df <- data.frame(year = 1:ncol(vel_discr), vel = vel_discr[j, ])
      vel_lm <- lm(vel ~ year, data = na.omit(vel_discr_df))
      vel_pred <- predict(vel_lm, newdata = data.frame(year = years))
      vel_discr_fin[j] <- vel_pred
    }

  }

  if (!correct_velocity_discrepancy) {
    vel_discr_mat <- matrix(0, nrow = length(domain), ncol = years)
  }

  ## Plot the discrepancy each year
  png(file = paste0(plot_dir, "pred_discrepancy_", data_date, ".png"), width = 1000, height = 1500, res = 150)

  par(mfrow = c(2, 1))
  matplot(domain / 1000, se_discr_mat,
    type = "l", lty = 1, col = rainbow(years - 1),
    main = "Surface elevation model discrepancy by year",
    xlab = "Distance along flowline (km)", ylab = "Discrepancy (m)"
  )
  # lines(domain/1000, se_discr_fin, col = "black", lwd = 2)

  matplot(domain / 1000, vel_discr_mat,
    type = "l", lty = 1, col = rainbow(years - 1),
    main = "Velocity model discrepancy by year",
    xlab = "Distance along flowline (km)", ylab = "Discrepancy (m/yr)"
  )
  # lines(domain/1000, vel_discr_fin, col = "black", lwd = 2)
  dev.off()

  ## Add discrepancy to simulated observations
  for (s in 1:length(prior_pred$results)) {
    prior_pred_obs[s, , , 1] <- prior_pred_obs[s, , , 1] + se_discr_mat
    prior_pred_obs[s, , , 2] <- prior_pred_obs[s, , , 2] + vel_discr_mat
  }

  for (s in 1:length(post_pred$results)) {
    # s <- 1
    post_pred_obs[s, , , 1] <- post_pred_obs[s, , , 1] + se_discr_mat
    post_pred_obs[s, , , 2] <- post_pred_obs[s, , , 2] + vel_discr_mat
  }
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
  domain = rep(domain / 1000, years),
  year = rep(2010 + 0:(years - 1), each = J),
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
  domain = rep(domain / 1000, years),
  year = rep(2010 + 0:(years - 1), each = J),
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

# gl_df <- data.frame(x = gl_obs/1000)

png(filename = paste0(plot_dir, "post_pred_se_by_year.png"), width = 2000, height = 3000, res = 200)
p1 <- ggplot(se_post_pred_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.4) +
  # geom_line(aes(y = prior_mean), color = "blue", size = 1) +
  geom_line(aes(y = obs), color = "black", lwd = 0.8) +
  # geom_vline(data = gl_df, aes(xintercept = x), lty = 2) + # plot GL position
  facet_wrap(~year, ncol = 2) +
  xlim(0, gl_obs / 1000) +
  # ylim(150, 1250) +
  labs(
    # title = "Prior vs posterior predictive surface elevation by year",
    y = "Surface elevation (m)", x = "Distance along flowline (km)"
  ) +
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
    labs(title = paste0("Year ", yr), y = "Surface elevation (m)", x = "Distance along flowline (km)") +
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
  facet_wrap(~year, ncol = 2, scales = "free_y") +
  xlim(0, 150) +
  # ylim(0, 4000) +
  labs(
    # title = "Prior vs posterior predictive surface velocity",
    y = "Surface velocity (m/yr)", x = "Distance along flowline (km)"
  ) +
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
    labs(title = paste0("Year ", yr), y = "Velocity (m/yr)", x = "Distance along flowline (km)") +
    theme_bw() +
    theme(text = element_text(size = 20))

  print(p_yr)

  dev.off()
}

## Compute RMSE
# rmse <- function(obs, sim) {
#   sqrt(mean((obs - sim)^2, na.rm = T))
# }
prior_se_rmse <- prior_vel_rmse <- c()
post_se_rmse <- post_vel_rmse <- c()

for (yr in 1:years) {
  ### Extract all final-year surface elevation obs from posterior predictive simulations
  obs_se_final <- surf_elev_data[, yr]
  obs_vel_final <- velocity_data[, yr]

  ### Discard observation locations with missing data
  valid_locs <- which(!is.na(obs_se_final))
  prior_pred_se_mean <- prior_pred_obs_mean[valid_locs, yr, 1]
  post_pred_se_mean <- post_pred_obs_mean[valid_locs, yr, 1]
  obs_se_final <- na.omit(obs_se_final)

  prior_pred_vel_mean <- prior_pred_obs_mean[valid_locs, yr, 2]
  post_pred_vel_mean <- post_pred_obs_mean[valid_locs, yr, 2]
  obs_vel_final <- obs_vel_final[valid_locs]

  prior_se_rmse[yr] <- rmse(obs_se_final, prior_pred_se_mean)
  post_se_rmse[yr] <- rmse(obs_se_final, post_pred_se_mean)
  prior_vel_rmse[yr] <- rmse(obs_vel_final, prior_pred_vel_mean)
  post_vel_rmse[yr] <- rmse(obs_vel_final, post_pred_vel_mean)

}

## CRPS

prior_crps_se_fin <- prior_crps_vel_fin <- c()
post_crps_se_fin <- post_crps_vel_fin <- c()
for (yr in 1:years) {
  ### Extract all final-year surface elevation obs from posterior predictive simulations
  prior_pred_se_final <- prior_pred_obs[, , yr, 1]
  post_pred_se_final <- post_pred_obs[, , yr, 1]
  obs_se_final <- surf_elev_data[, yr]

  ## Velocity
  prior_pred_vel_final <- prior_pred_obs[, , yr, 2]
  post_pred_vel_final <- post_pred_obs[, , yr, 2]
  obs_vel_final <- velocity_data[, yr]

  ### Discard observation locations with missing data
  valid_locs <- which(!is.na(obs_se_final))
  
  prior_pred_se_final <- prior_pred_se_final[, valid_locs]
  post_pred_se_final <- post_pred_se_final[, valid_locs]
  obs_se_final <- obs_se_final[valid_locs]

  prior_pred_vel_final <- prior_pred_vel_final[, valid_locs]
  post_pred_vel_final <- post_pred_vel_final[, valid_locs]
  obs_vel_final <- obs_vel_final[valid_locs]

  prior_crps_se_fin[yr] <- crps_multivariate(
    pred_samples = prior_pred_se_final,
    observed = obs_se_final
  )
  post_crps_se_fin[yr] <- crps_multivariate(
    pred_samples = post_pred_se_final,
    observed = obs_se_final
  )

  prior_crps_vel_fin[yr] <- crps_multivariate(
    prior_pred_vel_final,
    obs_vel_final
  )
  post_crps_vel_fin[yr] <- crps_multivariate(
    post_pred_vel_final,
    obs_vel_final
  )


}


## Collect RMSE into a data frame
rmse_df <- data.frame(
  year = as.integer(2010 + 0:(years - 1)),
  prior_se_RMSE = prior_se_rmse,
  post_se_RMSE = post_se_rmse,
  prior_vel_RMSE = prior_vel_rmse,
  post_vel_RMSE = post_vel_rmse
)

rmse_avg <- data.frame(
  year = c("Average"),
  prior_se_RMSE = mean(prior_se_rmse),
  post_se_RMSE = mean(post_se_rmse),
  prior_vel_RMSE = mean(prior_vel_rmse, na.rm = T),
  post_vel_RMSE = mean(post_vel_rmse, na.rm = T)
)

rmse_df_all <- rbind(rmse_df, rmse_avg)

## Save as a .csv file
write.csv(rmse_df_all, file = paste0(plot_dir, "rmse_summary_", data_date, ".csv"), row.names = F)

## Collect CRPS into a data frame
crps_df <- data.frame(
  year = as.integer(2010 + 0:(years - 1)),
  prior_se_CRPS = prior_crps_se_fin,
  post_se_CRPS = post_crps_se_fin,
  prior_vel_CRPS = prior_crps_vel_fin,
  post_vel_CRPS = post_crps_vel_fin
)

crps_avg <- data.frame(
  year = c("Average"),
  prior_se_CRPS = mean(prior_crps_se_fin),
  post_se_CRPS = mean(post_crps_se_fin),
  prior_vel_CRPS = mean(prior_crps_vel_fin, na.rm = T),
  post_vel_CRPS = mean(post_crps_vel_fin, na.rm = T)
)
crps_df_all <- rbind(crps_df, crps_avg)

## Save as a .csv file
write.csv(crps_df_all, file = paste0(plot_dir, "crps_summary_", data_date, ".csv"), row.names = F)

## Plot RMSE over time
png(filename = paste0(plot_dir, "rmse_plot_", data_date, ".png"), width = 1000, height = 800, res = 150)

se_rmse_plot <- ggplot() + 
  geom_line(data = rmse_df, aes(x = year, y = prior_se_RMSE), col = "blue") +
  geom_line(data = rmse_df, aes(x = year, y = post_se_RMSE), col = "red") +
  labs(title = "Surface elevation RMSE over time", y = "RMSE", x = "Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

vel_rmse_plot <- ggplot() + 
  geom_line(data = rmse_df, aes(x = year, y = prior_vel_RMSE), col = "blue") +
  geom_line(data = rmse_df, aes(x = year, y = post_vel_RMSE), col = "red") +
  labs(title = "Velocity RMSE over time", y = "RMSE", x = "Year") +
  theme_bw() +
  theme(text = element_text(size = 18))

grid.arrange(grobs = list(se_rmse_plot, vel_rmse_plot), nrow = 2, ncol = 1)
dev.off()

## Plot CRPS over time
png(filename = paste0(plot_dir, "crps_plot_", data_date, ".png"), width = 1000, height = 800, res = 150)  
se_crps_plot <- ggplot() + 
  geom_line(data = crps_df, aes(x = year, y = prior_se_CRPS), col = "blue") +
  geom_line(data = crps_df, aes(x = year, y = post_se_CRPS), col = "red") +
  labs(title = "Surface elevation CRPS over time", y = "CRPS", x = "Year") +
  theme_bw() +
  theme(text = element_text(size = 18))
vel_crps_plot <- ggplot() + 
  geom_line(data = crps_df, aes(x = year, y = prior_vel_CRPS), col = "blue") +
  geom_line(data = crps_df, aes(x = year, y = post_vel_CRPS), col = "red") +
  labs(title = "Velocity CRPS over time", y = "CRPS", x = "Year") +
  theme_bw() +
  theme(text = element_text(size = 18))
grid.arrange(grobs = list(se_crps_plot, vel_crps_plot), nrow = 2, ncol = 1)
dev.off()


## Plot prior vs posterior mean
fin_df <- se_post_pred_df %>% filter(year == 2020, domain %in% domain[valid_locs])
fin_df_vel <- vel_post_pred_df %>% filter(year == 2020, domain %in% domain[valid_locs])
png(filename = paste0(plot_dir, "prior_vs_posterior_mean_", data_date, ".png"), width = 1000, height = 1500, res = 150)
par(mfrow = c(2, 1))
plot(fin_df$obs, type = "l")
lines(fin_df$prior_mean, col = "blue")
lines(fin_df$post_mean, col = "red")
legend("bottomleft", legend = c("Observations", "Prior mean", "Posterior mean"), col = c("black", "blue", "red"), lty = 1)

plot(fin_df_vel$obs, type = "l")
lines(fin_df_vel$prior_mean, col = "blue")
lines(fin_df_vel$post_mean, col = "red")
legend("topleft", legend = c("Observations", "Prior mean", "Posterior mean"), col = c("black", "blue", "red"), lty = 1)

dev.off()

## Plot prior vs posterior spread
png(file = paste0(plot_dir, "prior_vs_posterior_pred_se_final_", data_date, ".png"), width = 1000, height = 1500, res = 150)

par(mfrow = c(2, 1))
matplot(t(prior_pred_se_final),
  type = "l", col = "lightblue",
  main = "Prior vs posterior predictive surface elevation samples",
  ylab = "Surface elevation (m)", xlab = "Location index"
)
matlines(t(post_pred_se_final), col = "lightpink", alpha = 0.25)
lines(obs_se_final, col = "black", lwd = 2)
lines(colMeans(prior_pred_se_final), col = "blue", lwd = 2)
lines(colMeans(post_pred_se_final), col = "red", lwd = 2)
legend("topright", legend = c("Prior samples", "Posterior samples", "Observations"), col = c("lightblue", "salmon", "black"), lty = 1)

## Same for the velocity
matplot(t(prior_pred_vel_final),
  type = "l", col = "lightblue",
  main = "Prior vs posterior predictive surface velocity samples",
  ylab = "Surface velocity (m/yr)", xlab = "Location index"
)
matlines(t(post_pred_vel_final), col = "lightpink", alpha = 0.25)
lines(obs_vel_final, col = "black", lwd = 2)
lines(colMeans(prior_pred_vel_final), col = "blue", lwd = 2)
lines(colMeans(post_pred_vel_final), col = "red", lwd = 2)
legend("topright", legend = c("Prior samples", "Posterior samples", "Observations"), col = c("lightblue", "salmon", "black"), lty = 1)

dev.off()