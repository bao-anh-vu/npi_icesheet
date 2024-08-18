### Main file ###

setwd("/home/babv971/SSA_model/EnKF")

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
library("mvnfast")
library("splines")
# library("fda")

source("surface_elev.R")
source("create_params.R")
source("create_ref.R")
source("solve_ssa_nl.R")
source("solve_velocity_azm.R")
source("solve_thickness.R")
source("mvnorm_sample.R")
source("get_obs.R")
source("simulate_bed.R")
source("simulate_friction.R")
source("azm_cond_sim.R")
source("initialise_ens.R")
source("ssa_plot_ini_ens.R")
source("propagate.R")
source("obs_operator.R")
source("run_enkf2.R")
source("run_bg_ens_backup.R")
# source("run_pf2.R")
source("run_pf_state_aug.R")
source("construct_bed_basis.R")
source("initialise_ice_thickness.R")
source("compute_block_weights.R")
source("create_pf_taper.R")
source("create_pf_smoother.R")

# source("ssa_enkf_plots.R")

## Seed for generating bed
ssa_seed <- 123
set.seed(ssa_seed)

## Some flags
add_process_noise <- T
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
use_ref_velocity <- F # use reference velocity as the initial previous velocity
use_velocity_dependent_R <- FALSE
run_analysis <- T
regenerate_ini_ens <- T # generate a fresh initial ensemble
use_basis_functions <- T
use_reference_thickness <- F
use_reference_bed <- F
use_reference_friction <- F
run_EnKF <- T
run_PF <- F
run_bg <- T
use_cov_taper <- TRUE # use covariance taper
use_variable_block <- T
smooth_by_weights <- F

## Presets
data_date <- "20220406b" # "20220607" 
Ne <- 50 # Ensemble size
years <- 2
steps_per_yr <- 25
n_beds <- 20
n_bed_obs <- 50
smoothing_factor <- 0.25
transformation <- "log10"

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

# 1. Reference run
reference <- create_ref(use_stored_steady_state = T,
                        use_stored_reference = T,
                        add_process_noise_in_ref = add_process_noise, 
                        rewrite_reference = F,
                        rewrite_steady_state = F,
                        data_date = data_date
)

## Extract some components for convenience
domain <- reference$domain
J <- length(domain)

# 2. Generate observatiobns of the top surface elevation, horizontal velocity and bed elevation
observations <- get_obs(reference, data_date = data_date, n_obs = n_bed_obs, 
                        get_new_obs = F, rewrite_obs = F) 

############################# Covariance tapering ##############################
wendland <- function(theta, D) {
  R <- D / theta
  W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
}

taper_mat <- wendland(0.1 * domain[length(domain)], rdist(domain, domain))
ens_taper1 <- kronecker(matrix(rep(1, 6), 3, 2), taper_mat)
ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)

if (!(run_EnKF && run_PF && run_bg)) {
  ini_ens_list <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", data_date, ".rds", sep = ""))
}

################################################################################
##                             EnKF implementation                            ##
################################################################################
if (run_EnKF) {
  
  # 3. Initialise ensemble
  print("Initialising ensemble...")
  
  if (regenerate_ini_ens) {
    
    ## Bed ##
    if (use_reference_bed) {
      ini_beds <- matrix(rep(reference$bedrock, Ne), J, Ne)
      
    } else if (use_basis_functions) {
      em_output <- construct_bed_basis(domain, bed_obs = observations$bed_obs,
                                       obs_locations = observations$bed_obs_locations,
                                       ini_params = c(-500, 0.5, 1000), n_basis = 50,
                                       n_beds = n_beds,
                                       parameterise_K = TRUE, use_sim_study = FALSE)
      ini_beds <- em_output$bed_ens
      
    } else {
      ini_beds <- simulate_bed(nsim = Ne,
                               domain = domain,
                               obs_locations = observations$bed_obs_locations,
                               obs = observations$bed_obs)
    }
    
    ## Friction ##
    print("Simulating friction coefficient...")
    if (use_reference_friction) {
      ini_friction <- matrix(rep(log10(reference$friction_coef), Ne), J, Ne)
    } else {
      simulated_alpha <- log10(simulate_friction(nsim = Ne, domain = reference$domain))
      ini_friction <- ifelse(is.na(simulated_alpha), 6, simulated_alpha) # replace NA values with 6
    }
    
    ## Now put the ensemble together
    ini_ens_list <- list()
    for (i in 1:n_beds) { # create 20 ensembles from 20 beds
      
      ## Ice thickness ##
      if (use_reference_thickness) {
        ini_thickness <- matrix(rep(reference$ini_thickness, Ne), J, Ne)
      } else {
        ini_thickness <- initialise_ice_thickness(domain,
                                                  n_sims = Ne, 
                                                  surface_obs = observations$surface_obs[, 1], # use observed z at t = 0
                                                  bed = ini_beds[, i], 
                                                  condsim_shelf = F)
      }
      
      ini_ens_list[[i]] <- rbind(ini_thickness,
                                 matrix(rep(ini_beds[, i], Ne), J, Ne),
                                 ini_friction)
    }
    # ini_ens <- rbind(ini_thickness, matrix(rep(ini_beds[, 1], Ne), J, Ne), ini_friction)
    
  } else {
    ini_ens_list <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", data_date, ".rds", sep = ""))
    # ini_ens_list <- ini_ens_list[c(1:Ne)]
  }
  
  ## Try plotting one of the ensembles
  ini_ens <- ini_ens_list[[1]]
  plot_ini_ens(ini_ens, reference, observations, print_friction_plot = TRUE)
  ################################ Initial velocity ##############################
  
  ini_velocity <- NULL
  
  if (use_ref_velocity) { # use reference velocity as initial velocity
    ini_velocity <- matrix(rep(reference$ini_velocity, Ne), J, Ne) 
  } else { # use a smoothed version of the observed velocity as initial velocity
    velocity_df <- data.frame(x = domain, u = observations$velocity_obs[, 1])
    smoothed_velocity <- loess(u ~ x, data = velocity_df, span = 0.05)$fitted
    ini_velocity <- matrix(rep(smoothed_velocity, Ne), J, Ne)
  }
  
  ## Run EnKF
  enkf1 <- proc.time()
  
  enkf_ens <- list()
  fin_ens <- list()
  fin_velocities <- list()
  
  for (i in 1:length(ini_ens_list)) {
    enkf_out <- run_enkf(domain, years, steps_per_yr, ini_ens_list[[i]][, 1:Ne], 
                         ini_velocity, observations, 
                         run_analysis = run_analysis,
                         add_process_noise = add_process_noise,
                         transformation = transformation) 
    
    enkf_ens[[i]] <- enkf_out$ens
    fin_ens[[i]] <- enkf_out$ens[[years + 1]] #enkf_ens <- enkf_out$ens
    # enkf_velocity_means <- enkf_out$velocity_means
    fin_velocities[[i]] <- enkf_out$velocities[[years + 1]] # enkf_velocities <- enkf_out$velocities
    
    enkf2 <- proc.time()
  } 
} else {
  fin_ens <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_ens_", data_date, ".rds", sep = ""))
  fin_velocities <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_velocities_", data_date, ".rds", sep = ""))
}

################################################################################
##                             PF implementation                              ##
################################################################################
error_ind <- rep(1, length(ini_ens_list))

if (run_PF) {
  pf1 <- proc.time()
  
  fin_pf_ens <- list()
  fin_pf_velocities <- list()
  fin_pf_weights <- list()
  
  for (i in 1:length(ini_ens_list)) {
    cat("PF bed number", i, "\n")
    # skip_to_next <- FALSE
    
    pf_out <- NULL
    attempt <- 0
    while (is.null(pf_out) && attempt < 2) {
      attempt <- attempt + 1
      cat("Attempt", attempt, "\n")
      tryCatch( expr = {
        pf_out <- run_pf(domain, years, steps_per_yr,
                         ini_ens_list[[i]], ini_velocity,
                         observations, run_analysis = run_analysis,
                         add_process_noise = add_process_noise, 
                         use_variable_block = use_variable_block,
                         lengthB_ground = 5,
                         lengthB_float = 3,
                         parallelise = T,
                         smooth_by_weights = smooth_by_weights)
        
        if (!is.null(pf_out)) {
          fin_pf_ens[[i]] <- pf_out$particles[[years + 1]] #enkf_ens <- enkf_out$ens
          fin_pf_weights[[i]] <- pf_out$weights[[years + 1]] #enkf_ens <- enkf_out$ens
          fin_pf_velocities[[i]] <- pf_out$velocities[[years + 1]] # enkf_velocities <- enkf_out$velocities
        }
        
        error_ind[i] <- 0
      },

      error = function(e) { cat("Error at bed number", i, "\n") }

      )

    }
    
    pf2 <- proc.time()
  }
  
} else {
  fin_pf_ens <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_ens_", data_date, ".rds", sep = ""))
  fin_pf_velocities <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_velocities_", data_date, ".rds", sep = ""))
}

################################################################################
##                              Background ensemble                           ##
################################################################################

fin_bg_ens <- list()
fin_bg_velocities <- list()

if (run_bg) {
  bg1 <- proc.time()
  
  for (i in 1:length(ini_ens_list)) { # parallelise later
    bg_out <- run_bg_ens(domain, years, steps_per_yr = 26, ini_ens_list[[i]], ini_velocity,
                         add_process_noise = add_process_noise)
    fin_bg_ens[[i]] <- bg_out$bg_ens[[years + 1]]
    fin_bg_velocities[[i]] <- bg_out$bg_velocities[[years + 1]]
  } 
  
  bg2 <- proc.time()
  
} else {
  fin_bg_ens <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_ens_", data_date, ".rds", sep = ""))
  fin_bg_velocities <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_velocities_", data_date, ".rds", sep = ""))
}

################################################################################
##                      Plotting the combined ensemble                        ##
################################################################################

par(mfrow = c(1,1))
combined_enkf_ens <- do.call(cbind, fin_ens)
combined_pf_ens <- do.call(cbind, fin_pf_ens)
combined_bg_ens <- do.call(cbind, fin_bg_ens)

combined_enkf_velocities <- do.call(cbind, fin_velocities)
combined_pf_velocities <- do.call(cbind, fin_pf_velocities)
combined_bg_velocities <- do.call(cbind, fin_bg_velocities)

## Grounding line ##
scale_factor <- 1 / J * domain[length(domain)] / 1000 # km scale

# Reference grounding line
ref.GL <- reference$grounding_line[(years) * steps_per_yr]

# EnKF grounding line
enkf.GL <- gl_migrate(H = rowMeans(combined_enkf_ens[1:J, ]), 
                      b = rowMeans(combined_enkf_ens[(J+1):(2*J), ])) * scale_factor

# PF grounding line
# pf.GL <- gl_migrate(H = rowMeans(combined_pf_ens[1:J, ]), 
                    # b = rowMeans(ini_ens_list[[1]])) * scale_factor

# Background GL
bg.GL <- gl_migrate(H = rowMeans(combined_bg_ens[1:J, ]), 
                    b = rowMeans(combined_bg_ens[(J+1):(2*J), ])) * scale_factor


## Ice thickness plot
if (plot_ice_thickness) {
  matplot(domain/1000, combined_bg_ens[1:J, ], type = "l", col = "lightgrey", 
          xlab = "Domain (km)", ylab = "Ice thickness (m)")
  matlines(domain/1000, combined_enkf_ens[1:J, ], lty = 2, col = "lightpink")
  # matlines(domain/1000, combined_pf_ens[1:J, ], lty = 2, col = "skyblue2")
  # lines(domain/1000, bg_state_means[1:J, years], col = "grey")
  lines(domain/1000, reference$all_thicknesses[, years+1], col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(combined_enkf_ens[1:J, ]), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_ens[1:J, ]), col = "royalblue", lwd = 1.5)
  abline(v = bg.GL, col = "grey", lty = 1)
  abline(v = enkf.GL, col = "red", lty = 2)
  # abline(v = pf.GL, col = "royalblue", lty = 3)
  abline(v = ref.GL, col = "black", lty = 2)
  legend("topright", legend = c("True thickness", "Background ensemble", 
                                "EnKF ensemble", "EnKF mean", "PF ensemble", "PF mean"),
         col = c("black", "grey", "pink", "red", "skyblue2", "royalblue"), 
         lty = 1, cex = 0.7)
}

## Velocity plot
if (plot_velocity) {
  matplot(domain/1000, combined_bg_velocities[1:J, ], type = "l", col = "lightgrey", 
          xlab = "Domain (km)", ylab = "Velocity (m/a)")
  # matlines(domain/1000, combined_pf_velocities[1:J, ], lty = 2, col = "skyblue2")
  matlines(domain/1000, combined_enkf_velocities[1:J, ], lty = 2, col = "lightpink")
  # lines(domain/1000, bg_state_means[1:J, years], col = "grey")
  lines(domain/1000, reference$all_velocities[, years+1], col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(combined_enkf_velocities[1:J, ]), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_velocities[1:J, ]), col = "royalblue", lwd = 1.5)
  abline(v = bg.GL, col = "grey", lty = 1)
  abline(v = enkf.GL, col = "red", lty = 2)
  # abline(v = pf.GL, col = "royalblue", lty = 3)
  abline(v = ref.GL, col = "black", lty = 2)
  legend("bottomright", legend = c("True thickness", "Background ensemble", 
                                   "EnKF ensemble", "EnKF mean", "PF ensemble", "PF mean"),
         col = c("black", "grey", "pink", "red", "skyblue2", "royalblue"), 
         lty = 1, cex = 0.7)
}

## Bed plot
if (plot_bed) {
  matplot(domain/1000, combined_bg_ens[(J+1):(2*J), ], type = "l", col = "lightgrey", 
          xlab = "Domain (km)", ylab = "Bed elevation (m)")
  matlines(domain/1000, combined_enkf_ens[(J+1):(2*J), ], col = "lightpink")
  # matlines(domain/1000, combined_pf_ens[(J+1):(2*J), ], col = "skyblue2")
  lines(domain/1000, reference$bedrock, col = "black", lwd = 1.5)
  lines(domain/1000, rowMeans(combined_enkf_ens[(J+1):(2*J), ]), col = "red", lwd = 1.5)
  # lines(domain/1000, rowMeans(combined_pf_ens[(J+1):(2*J), ]), col = "royalblue", lwd = 1.5)
  
  abline(v = bg.GL, col = "grey", lty = 1)
  abline(v = enkf.GL, col = "red", lty = 2)
  # abline(v = pf.GL, col = "royalblue", lty = 3)
  abline(v = ref.GL, col = "black", lty = 2)
  legend("bottomleft", legend = c("True bed", #"Background ensemble", 
                                  "EnKF ensemble", "EnKF mean"), 
                                  # "PF ensemble", "PF mean"),
         col = c("black", "pink", "red", "skyblue2", "royalblue"), 
         lty = 1, cex = 0.7)
}

## Recontruct covariance matrix
# fin_covmat <- 1/(Ne-1) * tcrossprod(10^combined_enkf_ens[1:J + 2*J])
# fin_mean_friction <- rowMeans(10^combined_enkf_ens[1:J + 2*J, ])
# lower_ci <- fin_mean_friction + qnorm(0.025) * sqrt(diag(fin_covmat))
# upper_ci <- fin_mean_friction + qnorm(0.975) * sqrt(diag(fin_covmat))

lower_ci <- c()
upper_ci <- c()

## Pointwise 95% CI
for (j in 1:J) {
  sorted <- sort(10^combined_enkf_ens[2*J + j, ])
  # take the 25th and 975th entries
  index_25 <- floor(ncol(combined_enkf_ens) / 100 * 2.5)
  index_975 <- ceiling(ncol(combined_enkf_ens) / 100 * 97.5)
  lower_ci[j] <- sorted[index_25]
  upper_ci[j] <- sorted[index_975]
}

# Friction plot
if (plot_friction) {
  plot_region <- 200:400#1:(ref.GL/domain[length(domain)] * 1000 * J)#(round(enkf.GL) - floor(J/4)):(round(enkf.GL) + 100)#400:950 #(2*J+1):(3*J)
  # matplot(domain[plot_region]/1000, 10^combined_bg_ens[2*J+plot_region, ], type = "l", col = "lightgrey",
          # xlab = "Domain (km)", ylab = "Friction coefficient (unit)")
  # matlines(domain[plot_region]/1000, 10^combined_enkf_ens[2*J+plot_region, ], col = "lightpink")
  # matlines(domain[plot_region]/1000, 10^combined_pf_ens[2*J+plot_region, ], col = "skyblue2")
  plot(domain[plot_region]/1000, rowMeans(10^combined_enkf_ens[2*J+plot_region, ]), 
       type = "l", col = "red", lwd = 1.5, ylab = "Friction coef (unit)", xlab = "Domain (km)")
  lines(domain[plot_region]/1000, reference$friction_coef[plot_region], col = "black", lwd = 1.5)
  # lines(domain[plot_region]/1000, rowMeans(fin_mean_friction[plot_region]), col = "red", lwd = 1.5)
  # lines(domain[plot_region]/1000, fin_mean_friction[plot_region], col = "red", lwd = 1.5)
  lines(domain[plot_region]/1000, lower_ci[plot_region], col = "salmon", lty = 2)
  lines(domain[plot_region]/1000, upper_ci[plot_region], col = "salmon", lty = 2)
  
  # lines(domain[plot_region]/1000, rowMeans(10^combined_pf_ens[2*J+plot_region, ]), col = "royalblue", lwd = 1.5)
  abline(v = bg.GL, col = "grey", lty = 1)
  abline(v = enkf.GL, col = "red", lty = 2)
  # abline(v = pf.GL, col = "royalblue", lty = 3)
  abline(v = ref.GL, col = "black", lty = 2)
  legend("topright", legend = c("True friction coefficient", "EnKF mean", "95% interval"),
         col = c("black", "red", "salmon"), lty = c(1,1,2), cex = 0.7)
}

## RMSE (time series)
rmse <- function(estimated, true) {
  stopifnot(length(estimated) == length(true))
  sum(sqrt((estimated - true)^2))
}

# if (!(run_EnKF && run_PF && run_bg)) { 
#   rmse.df <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/rmse_df_", data_date, ".rds", sep = ""))
# } else {
#   
#   enkf.rmse <- rmse(rowMeans(combined_enkf_ens[1:J, ]), reference$all_thicknesses[, years+1])
#   # pf.rmse <- rmse(rowMeans(combined_pf_ens[1:J, ]), reference$all_thicknesses[, years+1])
#   bg.rmse <- rmse(rowMeans(combined_bg_ens[1:J, ]), reference$all_thicknesses[, years+1])
#   
#   
# }

if (run_EnKF) {
  enkf.rmse <- rmse(rowMeans(combined_enkf_ens[1:J, ]), reference$all_thicknesses[, years+1])
} else {
  enkf.rmse <- NULL
}

if (run_PF) {
  pf.rmse <- rmse(rowMeans(combined_pf_ens[1:J, ]), reference$all_thicknesses[, years+1])
} else {
  pf.rmse <- NULL
}

if (run_bg) {
  bg.rmse <- rmse(rowMeans(combined_bg_ens[1:J, ]), reference$all_thicknesses[, years+1])
} else {
  bg.rmse <- NULL
}

# rmse.df <- data.frame(label = c("EnKF", "PF", "Background"), 
                      # rmse = c(enkf.rmse, pf.rmse, bg.rmse))

# print(rmse.df)

## PF error list
if (!run_PF) {
  error_ind <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/pf_error_ind_", data_date, ".rds", sep = ""))
  cat("PF error:", error_ind)
} else {
  cat("PF error:", error_ind)
}

################################ Save output ###################################
# saveRDS(reference, file = paste("/home/babv971/SSA_model/EnKF/Output/reference_", data_date, ".rds", sep = ""))
# saveRDS(observations, file = paste("/home/babv971/SSA_model/EnKF/Output/observations_", data_date, ".rds", sep = ""))
# saveRDS(ini_ens_list, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_list_", data_date, ".rds", sep = ""))
# saveRDS(fin_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_enkf_velocities_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_velocities_", data_date, ".rds", sep = ""))
# saveRDS(fin_pf_weights, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_pf_weights_", data_date, ".rds", sep = ""))
# saveRDS(fin_bg_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_ens_", data_date, ".rds", sep = ""))
# saveRDS(fin_bg_velocities, file = paste("/home/babv971/SSA_model/EnKF/Output/fin_bg_velocities_", data_date, ".rds", sep = ""))
# saveRDS(error_ind, file = paste("/home/babv971/SSA_model/EnKF/Output/pf_error_ind_", data_date, ".rds", sep = ""))
# saveRDS(rmse.df, file = paste("/home/babv971/SSA_model/EnKF/Output/rmse_df_", data_date, ".rds", sep = ""))
