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
sets <- 51:100 # 6:20
# use_missing_pattern <- T
use_basal_melt_data <- F
correct_model_discrepancy <- T

setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

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

## Flowline data
# flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
flowline <- qread(paste0("./data/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

## Load physical parameters
params <- qread(file = paste0("./data/training_data/", "/phys_params_", data_date, ".qs"))

if (!use_basal_melt_data) {
    params$ab <- 0
}

## Posterior samples
post_fric_samples <- qread(file = paste0(pred_output_dir, "fric_samples_real_", data_date, ".qs"))
post_bed_samples <- qread(file = paste0(pred_output_dir, "bed_samples_real_", data_date, ".qs"))

## Initialise surface velocity at the 10th year
n_samples <- 2
se_grounded <- na.omit(surf_elev_data[, 10]) # Use surface elevation in the year 2000 to initialise ice thickness
H_ini <- se_grounded - post_bed_samples[1:length(se_grounded), 1:n_samples]
length_shelf <- J - length(se_grounded)

# H_gl <- - bed_sim[(gl_ind+1)] * params$rho_w / params$rho_i #- 100 # minus an offset to satisfy flotation condition 
H_shelf <- rep(500, length_shelf) #seq(from = H_gl, to = 500, length.out = length_shelf)
H_shelf_mat <- matrix(rep(H_shelf, n_samples), ncol = n_samples)
# thickness_at_gl <- - bed_sim[gl_ind][1] * params$rho_w / params$rho_i
# H_shelf <- seq(thickness_at_gl - 1, 500, length.out = length_shelf)

ini_thickness <- rbind(H_ini, H_shelf_mat)

z <- lapply(1:n_samples, function(s) {
  get_surface_elev(H = ini_thickness[, s], b = post_bed_samples[, s], 
  z0 = 0, rho = params$rho_i, rho_w = params$rho_w, include_GL = TRUE)
})

ini_velocity <- velocity_data[, 10] #+ rnorm(J, 0, 50)
# Take last non-NA value of velocity data
v_tail <- tail(ini_velocity[which(!is.na(ini_velocity))], 1)
ini_velocity[is.na(ini_velocity)] <- v_tail

ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))


post_param_list <- lapply(1:n_samples, function(r) {
  list(
    friction = post_fric_samples[, r], #* 1e6 * params$secpera^(1 / params$n),
    bedrock = post_bed_samples[, r], #+ bed_mean
#     ini_thickness = ssa_steady$current_thickness,
#   ini_velocity = ssa_steady$current_velocity
    ini_thickness = ini_thickness[, r],
    ini_velocity = ini_velocity #+ rnorm(J, 0,
  )
})

one_year_pred <- sim_obs(
  param_list = post_param_list,
  domain = flowline_dist,
  phys_params = params,
  years = 1, # sim_beds = T,
  smb = smb_avg,
  basal_melt = avg_melt_rate
  # log_transform = log_transform
)

## Plot posterior sample of bed and friction
png(paste0("./plots/temp/post_params_", data_date, ".png"), width = 1000, height = 1000, res = 100)
par(mfrow = c(2, 1))
plot(post_fric_samples[, 1], type = 'l', ylim = range(post_fric_samples), xlab = "Distance along flowline (m)", ylab = "Friction coefficient", main = "Posterior samples of friction coefficient")
plot(post_bed_samples[, 1], type = 'l', ylim = range(post_bed_samples), xlab = "Distance along flowline (m)", ylab = "Bed elevation (m)", main = "Posterior samples of bed elevation")
dev.off()
