### Main file ###

setwd("/home/babv971/SSA_model/CNN/simbed/")

rm(list = ls())

library("parallel")
library("Matrix")
library("qlcMatrix")
library("fastmatrix")
library("expm")
library("R.utils")
library("sp")
library("fields")
library("tidyr")
library("dplyr")
library("ggplot2")
library("matrixStats")
library("mvtnorm")
# library("mvnfast")
library("splines")
# library("fda")

# source("./source/enkf/create_params.R")
# source("./source/enkf/create_ref.R")
source("./source/solve_ssa_nl.R")
# source("./source/enkf/solve_velocity_azm.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
# source("./source/enkf/mvnorm_sample.R")
source("./source/enkf/get_obs.R")
# source("./source/enkf/simulate_bed.R")
# source("./source/enkf/simulate_friction.R")
# source("./source/enkf/azm_cond_sim.R")
source("./source/enkf/initialise_ens.R")
source("./source/enkf/ssa_plot_ini_ens.R")
source("./source/enkf/propagate.R")
source("./source/enkf/obs_operator.R")
source("./source/enkf/run_enkf.R")
source("./source/enkf/surface_elev.R")

# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
source("./source/enkf/initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# source("ssa_enkf_plots.R")

## Seed for generating bed
ssa_seed <- 123
set.seed(ssa_seed)

## Some flags
add_process_noise <- T
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
use_true_velocity <- T # use reference velocity as the initial previous velocity
# use_velocity_dependent_R <- FALSE
# run_analysis <- T
regenerate_ini_ens <- T # generate a fresh initial ensemble
use_basis_functions <- F
use_true_thickness <- T
use_true_bed <- T
use_true_friction <- T
run_EnKF <- T
# run_PF <- F
# run_bg <- F
use_cov_taper <- TRUE # use covariance taper
# use_variable_block <- T
# smooth_by_weights <- F

save_enkf_output <- T
save_bg_output <- F

## Presets
data_date <- "20220329" #"20230518"
output_date <- "20240320" #"20240518"
Ne <- 50 # Ensemble size
years <- 2#0 #40
steps_per_yr <- 52 #100
n_params <- 1 #20 #number of beds
# n_bed_obs <- 100
# smoothing_factor <- 0.5 # for the PF

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T 

# 1. Reference run
# print("Simulating ground truth...")
# reference <- create_ref(use_stored_steady_state = T,
#                         use_stored_reference = T,
#                         add_process_noise_in_ref = T,
#                         rewrite_steady_state = F, 
#                         rewrite_reference = F,
#                         data_date = data_date 
#                         )

# ## Extract some components for convenience
# domain <- reference$domain
# J <- length(domain)

# # 2. Generate observations of the top surface elevation, horizontal velocity and bed elevation
# print("Generating observations...")
# observations <- get_obs(reference, data_date = data_date, n_obs = n_bed_obs, 
#                         get_new_obs = T, rewrite_obs = F) 

## SSA model info
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
reference <- readRDS(file = paste("./training_data/initial_conds/reference_20220329.rds", sep = ""))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed and friction from NN output
sets <- 1:5 #10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./training_data/", setsf)
output_dir <- paste0("./output/all/", setsf)
test_data <- readRDS(file = paste0(data_dir, "/test_data_", output_date, ".rds"))

true_surface_elevs <- test_data$true_surface_elevs_test
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- test_data$true_fric

pred_fric <- readRDS(file = paste0(output_dir, "/pred_fric_", output_date, ".rds"))
pred_bed <- readRDS(file = paste0(output_dir, "/pred_bed_", output_date, ".rds"))
# pred_gl <- readRDS(pred_gl, file = paste0(output_dir, "/pred_gl_", output_date, ".rds"))

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1/3)
pred_fric <- pred_fric * fric_scale

## Read surface observations
# surface_obs <- test_data$input * test_data$input_sd + test_data$input_mean
surface_elev <- test_data$input[,,,1] * test_data$input_sd[1] + test_data$input_mean[1]
velocity <- test_data$input[,,,2] * test_data$input_sd[2] + test_data$input_mean[2]
# surface_elev <- surface_obs[,,,1]
# velocity <- surface_obs[,,,2]

################################
##      State inference       ##
################################

## Sample from test set and do state inference
s <- 1
fric_s <- pred_fric[, s]
bed_s <- pred_bed[, s]

# surface_obs_s <- surface_obs[s,,,]
surface_elev_s <- surface_elev[s,,]
velocity_s <- velocity[s,,]

surface_obs_list <- list(surface_elev = surface_elev_s, velocity = velocity_s)

# 3. Initialise ensemble
print("Initialising ensemble...")

if (regenerate_ini_ens) {
  
  ## Bed ##
  if (use_true_bed) {
    ini_beds <- matrix(rep(true_bed[s, ], Ne), J, Ne)
    
  } else {
    ini_beds <- matrix(rep(bed_s, Ne), J, Ne)
  }
  
  ## Friction ##
  print("Simulating friction coefficient...")
  if (use_true_friction) {
    ini_friction <- matrix(rep(log10(true_fric[s, ]*fric_scale), Ne), J, Ne)
  } else {
    ini_friction <- matrix(rep(log10(fric_s), Ne), J, Ne)
    # simulated_alpha <- log10(simulate_friction(nsim = Ne, domain = reference$domain))
    # ini_friction <- ifelse(is.na(simulated_alpha), 6, simulated_alpha) # replace NA values with 6
  }
  
  ## Now put the ensemble together
  print("Calculating ice thicknesses based on simulated beds and observed top surface elevation...")
  ini_ens_list <- list()
  for (i in 1:n_params) { # create 20 ensembles from 20 beds
    
    ## Ice thickness ##
    if (use_true_thickness) {
      ini_thickness <- matrix(rep(true_thicknesses[1,,1], Ne), J, Ne)
    } else {
      ini_thickness <- initialise_ice_thickness(domain, n_sims = Ne, 
                                                surface_obs = surface_elev_s[, 1], # use observed z at t = 0
                                                bed = ini_beds[, i], 
                                                condsim_shelf = F)
    }

    # ini_ens_list[[i]] <- rbind(ini_thickness,
    #                            matrix(rep(ini_beds[, i], Ne), J, Ne),
    #                            ini_friction)
    ini_ens <- rbind(ini_thickness,
                    matrix(rep(ini_beds[, i], Ne), J, Ne),
                    ini_friction)
  }
  # ini_ens <- rbind(ini_thickness, matrix(rep(ini_beds[, 1], Ne), J, Ne), ini_friction)
  
  if (save_enkf_output) {
    saveRDS(ini_ens_list, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", output_date, ".rds", sep = ""))
  }

} else {
  ini_ens_list <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", output_date, ".rds", sep = ""))
  # ini_ens <- ini_ens[, 1:Ne]
}

## Try plotting one of the ensembles
# param_ind <- 1
# ini_state_params <- rbind(ini_ens_list[[param_ind]], 
#                           ini_bed_list[[param_ind]], 
#                           ini_friction_list[[param_ind]])
# plot_ini_ens(ini_ens, reference, #observations, 
#               print_bed_plot = T, print_friction_plot = T)
################################ Initial velocity ##############################

ini_velocity <- NULL

if (use_true_velocity) { # use reference velocity as initial velocity
  ini_velocity <- matrix(rep(true_velocities[1,,1], Ne), J, Ne) 
} else { # use a smoothed version of the observed velocity as initial velocity
  velocity_df <- data.frame(x = domain, u = velocity_s[, 1])#u = observations$velocity_obs[, 1])
  smoothed_velocity <- loess(u ~ x, data = velocity_df, span = 0.05)$fitted
  ini_velocity <- matrix(rep(smoothed_velocity, Ne), J, Ne)
}

############################# Covariance tapering ##############################
# wendland <- function(theta, D) {
#   R <- D / theta
#   W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
# }
# 
# taper_mat <- wendland(0.1 * domain[length(domain)], rdist(domain, domain))
# ens_taper1 <- kronecker(matrix(rep(1, 6), 3, 2), taper_mat)
# ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)

################################################################################
##                             EnKF implementation                            ##
################################################################################
if (run_EnKF) {

    ## Process noise parameters
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 50e3
    R <- exp_cov(D, l)

    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    process_noise_info <- list(corrmat_chol = L, length_scale = l)


  enkf1 <- proc.time()
  
  # fin_enkf_ens <- list()
  # fin_enkf_beds <- list()
  # fin_enkf_friction <- list()
  # fin_velocities <- list()
  
  # for (i in 1:length(ini_ens_list)) { 
    # enkf_out <- run_enkf(domain, years, steps_per_yr, 
    #                      ini_ens = ini_ens_list[[i]][1:J, ],
    #                      ini_bed = ini_ens_list[[i]][J + 1:J, ],
    #                      ini_friction_coef = ini_ens_list[[i]][2*J + 1:J, ],
    #                      ini_velocity = ini_velocity, 
    #                      observations = surface_obs_list, 
    #                     #  observations = observations, 
    #                      run_analysis = run_analysis,
    #                      add_process_noise = add_process_noise)
    # enkf_means <- enkf_out$ens_means
    # fin_enkf_ens[[i]] <- enkf_out$ens#[[years + 1]] #enkf_ens <- enkf_out$ens
    # fin_enkf_beds[[i]] <- enkf_out$bed #enkf_ens <- enkf_out$ens
    # fin_enkf_friction[[i]] <- enkf_out$friction_coef #enkf_ens <- enkf_out$ens
    # # enkf_velocity_means <- enkf_out$velocity_means
    # fin_velocities[[i]] <- enkf_out$velocities#[[years + 1]] # enkf_velocities <- enkf_out$velocities
    
    enkf_out <- run_enkf(domain, years, steps_per_yr, 
                         ini_ens = ini_ens[1:J, ],
                         ini_bed = ini_ens[J + 1:J, ],
                         ini_friction_coef = ini_ens[2*J + 1:J, ],
                         ini_velocity = ini_velocity, 
                         observations = surface_obs_list, 
                        #  observations = observations, 
                        #  run_analysis = run_analysis,
                        #  add_process_noise = add_process_noise)
                         run_analysis = T,
                         add_process_noise = T,
                         process_noise_info = process_noise_info)
                         

    enkf_thickness <- enkf_out$ens
    enkf_bed <- enkf_out$bed
    enkf_friction <- enkf_out$friction_coef
    enkf_velocities <- enkf_out$velocities
    # enkf_means <- enkf_out$ens_means
    # fin_enkf_ens[[i]] <- enkf_out$ens#[[years + 1]] #enkf_ens <- enkf_out$ens
    # fin_enkf_beds[[i]] <- enkf_out$bed #enkf_ens <- enkf_out$ens
    # fin_enkf_friction[[i]] <- enkf_out$friction_coef #enkf_ens <- enkf_out$ens
    # # enkf_velocity_means <- enkf_out$velocity_means
    # fin_velocities[[i]] <- enkf_out$velocities#[[years + 1]] # enkf_velocities <- enkf_out$velocities
    
    enkf2 <- proc.time()
  # } 

  if (save_enkf_output) {
    # saveRDS(fin_enkf_ens, file = paste("./output/enkf/fin_enkf_ens_", output_date, ".rds", sep = ""))
    # saveRDS(fin_enkf_beds, file = paste("./output/enkf/fin_enkf_beds_", output_date, ".rds", sep = ""))
    # saveRDS(fin_enkf_friction, file = paste("./output/enkf/fin_enkf_friction_", output_date, ".rds", sep = ""))
    # saveRDS(fin_velocities, file = paste("./output/enkf/fin_enkf_velocities_", output_date, ".rds", sep = ""))

    saveRDS(enkf_thickness, file = paste0(output_dir, "/enkf_thickness_", output_date, ".rds", sep = ""))
    saveRDS(enkf_bed, file = paste0(output_dir, "/enkf_bed_", output_date, ".rds", sep = ""))
    saveRDS(enkf_friction, file = paste0(output_dir, "/enkf_friction_", output_date, ".rds", sep = ""))
    saveRDS(enkf_velocities, file = paste0(output_dir, "/enkf_velocities_", output_date, ".rds", sep = ""))
  }

} else {
    enkf_thickness <- saveRDS(enkf_thickness, file = paste0(output_dir, "/enkf_thickness_", output_date, ".rds", sep = ""))
    enkf_bed <- saveRDS(enkf_bed, file = paste0(output_dir, "/enkf_bed_", output_date, ".rds", sep = ""))
    enkf_friction <- saveRDS(enkf_friction, file = paste0(output_dir, "/enkf_friction_", output_date, ".rds", sep = ""))
    enkf_velocity <- saveRDS(enkf_velocities, file = paste0(output_dir, "/enkf_velocities_", output_date, ".rds", sep = ""))
}

################################################################################
##                      Plotting the combined ensemble                        ##
################################################################################

par(mfrow = c(1,1))
# combined_enkf_ens <- do.call(cbind, fin_enkf_ens)
# combined_enkf_beds <- do.call(cbind, fin_enkf_beds)
# combined_enkf_friction <- do.call(cbind, fin_enkf_friction)
# combined_enkf_velocities <- do.call(cbind, fin_velocities)


plot_dir <- paste0("./plots/all/", setsf)
## Ice thickness plot
if (plot_ice_thickness) {
  png(paste0(plot_dir, "/ice_thickness.png"), width = 2000, height = 1000, res = 300)
  # matplot(domain/1000, combined_bg_ens[1:J, ], type = "l", col = "lightgrey", 
  #         xlab = "Domain (km)", ylab = "Ice thickness (m)")
  matplot(domain/1000, enkf_thickness[[years+1]], lty = 2, type = "l", col = "lightpink")
  # matlines(domain/1000, combined_pf_ens[1:J, ], lty = 2, col = "skyblue2")
  # lines(domain/1000, rowMeans(combined_bg_ens[1:J, ]), col = "grey")
  # lines(domain/1000, reference$all_thicknesses[, years+1], col = "black", lwd = 1.5)
  lines(domain/1000, true_thicknesses[s, , years+1], col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(enkf_thickness[[years+1]]), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_ens[1:J, ]), col = "royalblue", lwd = 1.5)
  legend("topright", legend = c("True thickness", #"Background ensemble", 
                                "EnKF ensemble", "EnKF mean"), #, "PF ensemble", "PF mean"),
         col = c("black", #"grey", 
         "pink", "red"), #, "skyblue2", "royalblue"), 
         lty = 1, cex = 0.7)
  dev.off()
}

## Velocity plot
if (plot_velocity) {
  png(paste0(plot_dir, "/velocity.png"), width = 2000, height = 1000, res = 300)
  # matplot(domain/1000, combined_bg_velocities[1:J, ], type = "l", col = "lightgrey", 
  #         xlab = "Domain (km)", ylab = "Velocity (m/a)")
  # matlines(domain/1000, combined_pf_velocities[1:J, ], lty = 2, col = "skyblue2")
  matplot(domain/1000, enkf_velocities[[years+1]], lty = 2, type = "l", col = "lightpink")
  # lines(domain/1000, rowMeans(combined_bg_velocities[1:J, ]), col = "grey")
  # lines(domain/1000, reference$all_velocities[, years+1], col = "black", lwd = 1.5)
  lines(domain/1000, true_velocities[s, , years+1], col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(enkf_velocities[[years+1]]), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_velocities[1:J, ]), col = "royalblue", lwd = 1.5)
  legend("bottomright", legend = c("True thickness", 
                                    # "Background ensemble", 
                                   "EnKF ensemble", "EnKF mean"), #, "PF ensemble", "PF mean"),
         col = c("black", #"grey", 
         "lightpink", "red"), #, "skyblue2", "royalblue"), 
         lty = 1, cex = 0.7)
  dev.off()
}

## Bed plot
if (plot_bed) {
  png(paste0(plot_dir, "/bed.png"), width = 2000, height = 1000, res = 300)
  # matplot(domain/1000, combined_bg_beds, type = "l", col = "lightgrey", 
  #         xlab = "Domain (km)", ylab = "Elevation (m)")
  matplot(domain/1000, enkf_bed, type = "l", col = "lightpink",
          ylab = "Bed elevation (m)", xlab = "Domain (km)")
  # matlines(domain/1000, combined_pf_beds, col = "skyblue2")
  # lines(domain/1000, reference$bedrock, col = "black", lwd = 1.5)
  lines(domain/1000, true_bed[s, ], col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(enkf_bed), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_beds), col = "royalblue", lwd = 1.5)
  legend("bottomleft", legend = c("True bed", "Predicted bed"), #"Background ensemble", 
                                  # "EnKF ensemble", "EnKF mean"),
         col = c("black", #"grey", "pink", 
         "red"), lty = 1)
  dev.off()
}

# Friction plot
if (plot_friction) {
  png(paste0(plot_dir, "/friction.png"), width = 2000, height = 1000, res = 300)
  plot_region <- 1:J #500:1000 #(2*J+1):(3*J)
  # matplot(domain[plot_region]/1000, 10^combined_bg_friction[plot_region, ], type = "l", col = "lightgrey", 
  #         xlab = "Domain (km)", ylab = "Friction coefficient (unit)")
  matplot(domain[plot_region]/1000, 10^enkf_friction[plot_region, ]/fric_scale, 
          type = "l", col = "lightpink", 
          ylab = "Friction (unit)", xlab = "Domain (km)")
  # matlines(domain[plot_region]/1000, 10^combined_pf_friction[plot_region, ], col = "skyblue2")
  # lines(domain[plot_region]/1000, reference$friction_coef[plot_region], col = "black", lwd = 1.5)
  lines(domain/1000, true_fric[s, ], col = "black", lwd = 1.5)
  lines(domain[plot_region]/1000, rowMeans(10^enkf_friction[plot_region, ]/fric_scale), col = "red", lwd = 1.5)
  # lines(domain[plot_region]/1000, rowMeans(10^combined_pf_friction[plot_region, ]), col = "royalblue", lwd = 1.5)
  legend("bottomleft", legend = c("True friction", "Predicted friction"),
                                  # "Background ensemble", 
                                  # "PF ensemble", "PF mean", 
                                  # "EnKF ensemble", "EnKF mean"),
         col = c("black", #"lightgrey", 
                 # "skyblue2", "royalblue", 
                #  "lightpink", 
                 "red"), lty = 1, cex = 0.7)
  dev.off()
}

## RMSE (time series)
rmse <- function(estimated, true) {
  stopifnot(length(estimated) == length(true))
  sum(sqrt((estimated - true)^2))
}

enkf.rmse <- rmse(rowMeans(enkf_thickness[[years+1]]), true_thicknesses[s, , years+1])

# enkf.rmse <- rmse(rowMeans(combined_enkf_ens[1:J, ]), reference$all_thicknesses[, years+1])
# pf.rmse <- rmse(rowMeans(combined_pf_ens[1:J, ]), reference$all_thicknesses[, years+1])
# bg.rmse <- rmse(rowMeans(combined_bg_ens[1:J, ]), reference$all_thicknesses[, years+1])

# rmse.df <- data.frame(label = c("EnKF", "PF", "Background"), 
#                       rmse = c(enkf.rmse, pf.rmse, bg.rmse))

# rmse.df <- data.frame(label = c("EnKF", "Background"), 
#                       rmse = c(enkf.rmse, bg.rmse))
# print(rmse.df)

# ## PF error
# if (!run_PF) {
#   error_ind <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/pf_error_ind_", data_date, ".rds", sep = ""))
#   cat("PF error:", error_ind)
# }


################################ Save output ###################################
# saveRDS(reference, file = paste("/home/babv971/SSA_model/EnKF/Output/reference_", data_date, ".rds", sep = ""))
# saveRDS(observations, file = paste("/home/babv971/SSA_model/EnKF/Output/observations_", data_date, ".rds", sep = ""))
# saveRDS(ini_ens_list, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", data_date, ".rds", sep = ""))
# saveRDS(ini_bed_list, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_bed_list_", data_date, ".rds", sep = ""))
# saveRDS(ini_friction_list, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_friction_list_", data_date, ".rds", sep = ""))

# saveRDS(fin_enkf_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_enkf_beds, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_beds_", data_date, ".rds", sep = ""))
# saveRDS(fin_enkf_friction, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_friction_", data_date, ".rds", sep = ""))
# saveRDS(fin_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_velocities_", data_date, ".rds", sep = ""))
# 
# saveRDS(fin_pf_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_beds, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_beds_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_friction, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_friction_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_velocities_", data_date, ".rds", sep = ""))
# saveRDS(error_ind, file = paste("/home/babv971/SSA_model/EnKF/Output/pf_error_ind_", data_date, ".rds", sep = ""))
# 
# saveRDS(fin_bg_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_bg_beds, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_beds_", data_date, ".rds", sep = ""))
# saveRDS(fin_bg_friction, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_friction_", data_date, ".rds", sep = ""))
# saveRDS(fin_bg_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_velocities_", data_date, ".rds", sep = ""))
