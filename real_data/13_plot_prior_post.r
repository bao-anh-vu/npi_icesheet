## Plot prior vs posterior distributions of bed and friction

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

data_date <- "20241111" # "20241103"
sets <- 51:100 # 6:20
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
n_samples <- 100
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

prior_bed_mean <- apply(prior_bed_samples_basis, 1, mean)
prior_bed_lci <- apply(prior_bed_samples_basis, 1, quantile, probs = 0.025)
prior_bed_uci <- apply(prior_bed_samples_basis, 1, quantile, probs = 0.975)

prior_fric_mean <- apply(prior_fric_samples_basis, 1, mean)
prior_fric_lci <- apply(prior_fric_samples_basis, 1, quantile, probs = 0.025)
prior_fric_uci <- apply(prior_fric_samples_basis, 1, quantile, probs = 0.975)

post_bed_mean <- apply(post_bed_samples, 1, mean)
post_bed_lci <- apply(post_bed_samples, 1, quantile, probs = 0.025)
post_bed_uci <- apply(post_bed_samples, 1, quantile, probs = 0.975)

post_fric_mean <- apply(post_fric_samples, 1, mean)
post_fric_lci <- apply(post_fric_samples, 1, quantile, probs = 0.025)
post_fric_uci <- apply(post_fric_samples, 1, quantile, probs = 0.975)

bed_df <- data.frame(
  domain = domain/1000,
  prior_mean = as.vector(prior_bed_mean),
  prior_lci = as.vector(prior_bed_lci),
  prior_uci = as.vector(prior_bed_uci),
  post_mean = as.vector(post_bed_mean),
  post_lci = as.vector(post_bed_lci),
  post_uci = as.vector(post_bed_uci)
)

fric_df <- data.frame(
  domain = domain/1000,
  prior_mean = as.vector(prior_fric_mean),
  prior_lci = as.vector(prior_fric_lci),
  prior_uci = as.vector(prior_fric_uci),
  post_mean = as.vector(post_fric_mean),
  post_lci = as.vector(post_fric_lci),
  post_uci = as.vector(post_fric_uci)
)

## Plot prior vs posterior in ggplot
bed_p <- ggplot(bed_df, aes(x = domain)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
  geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.4) +
  geom_line(aes(y = prior_mean), color = "blue", lwd = 1) +
  geom_line(aes(y = post_mean), color = "red", lwd = 1) +
  xlim(0, 150) +
  labs(title = "Prior vs posterior distribution of bed", 
       y = "Bed elevation (m)", x = "Dist. along flowline (km)") +
  theme_minimal()

png(file = paste0(plot_dir, "prior_vs_post_bed_", data_date, ".png"), 
      width = 800, height = 600)
print(bed_p)
dev.off()   

## Do a generic plot of a function with uncertainty bands

x <- seq(0, 2*pi, length.out = 1000)
prior_mean <- 2 + 0.5 * sin(x)

prior_lci <- prior_mean - 0.4
prior_uci <- prior_mean + 0.4

post_mean <- 2 + 0.5 * sin(x + 0.5) - 0.05 * cos(1.5 * pi * x) + 0.08 * cos(1.2 * pi * x)
post_lci <- post_mean - 0.25
post_uci <- post_mean + 0.25

gen_df <- data.frame(
  x = x,
  prior_mean = prior_mean,
  prior_lci = prior_lci,
  prior_uci = prior_uci,
  post_mean = post_mean,
  post_lci = post_lci,
  post_uci = post_uci
)

gen_p <- ggplot(gen_df, aes(x = x)) +
  geom_ribbon(aes(ymin = prior_lci, ymax = prior_uci), fill = "lightblue", alpha = 0.75) +
#   geom_ribbon(aes(ymin = post_lci, ymax = post_uci), fill = "salmon", alpha = 0.4) +
  geom_line(aes(y = prior_mean), color = "blue", lwd = 1) +
#   geom_line(aes(y = post_mean), color = "red", lwd = 1) +
  labs(
    # title = "Prior vs posterior distribution", 
       y = "Friction", x = "Flowline (km)") +
  theme_bw()

png(file = paste0(plot_dir, "generic_prior_only", data_date, ".png"), 
      width = 800, height = 500, res = 150)
print(gen_p)
dev.off()


## Plot generic prior vs posterior

# Load required libraries
# library(ggplot2)
# library(dplyr)

# Define a sequence for theta
theta <- seq(0, 3, length.out = 1000)

# Define parameters
# a_prior <- prior_mean[1] #10
# b_prior <- post_mean[1] #10
# a_likelihood <- 15
# b_likelihood <- 7
# a_posterior <- a_prior + a_likelihood - 1
# b_posterior <- b_prior + b_likelihood - 1

# Create data frames
df <- data.frame(
  theta = theta,
  Prior = dnorm(theta, prior_mean[1], 0.2),
#   Likelihood = dbeta(theta, a_likelihood, b_likelihood),
#   Posterior = dnorm(theta, post_mean[1], 0.1)
)

# Reshape data for ggplot
df_long <- df %>%
  tidyr::pivot_longer(#cols = c(Prior, Likelihood, Posterior),
                    cols = c(Prior), #Posterior),
                      names_to = "Distribution",
                      values_to = "Density")

# Plot
dist_p <- ggplot(df_long, aes(x = theta, y = Density, fill = Distribution)) +
  geom_area(alpha = 0.6, position = 'identity') +
  scale_fill_manual(values = c("Prior" = "lightblue")) +                            #    "Likelihood" = "purple", 
                            #    "Posterior" = "salmon")) +
  geom_vline(xintercept = prior_mean[1], color = "blue", linetype = "dashed", size = 1) +
#   geom_vline(xintercept = post_mean[1], color = "red", linetype = "dashed", size = 1) +
  xlim(1.2, 2.8) +
  ylim(0, 4) +                           
  labs(x = "Friction", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 30))

png(file = paste0(plot_dir, "generic_bayes_update_prior_", data_date, ".png"), 
      width = 2000, height = 1000, res = 300)
print(dist_p)
dev.off()