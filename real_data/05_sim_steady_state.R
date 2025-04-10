setwd("~/SSA_model/CNN/real_data/")

library(dplyr)
library(ggplot2)
library(fields)
library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
library(R.utils)
# library("sp")
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(qs)
# library(sf)

source("./source/create_params.R")
source("./source/simulate_friction.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
# source("./source/fit_basis.R")

rerun_steady_state <- T
use_basal_melt_data <- T

data_dir <- "./data/"
data_date <- "20241111" #"20241103"

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
params$B <- 0.7 * 1e6 * params$secpera^params$m

## Flowline data
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
J <- nrow(flowline) # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## Simulate bed
# bed_sim <- create_bed(x = flowline_dist)
bed_sim <- qread(file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))

print("Simulating friction coefficient...")

## Simulate friction coefficient
# secpera <- 31556926 #seconds per annum
# n <- 3.0 # exponent in Glen's flow law
# m <- 1/n # friction exponent
L <- flowline_dist[J] - flowline_dist[1]
fric_sim <- create_fric_coef(flowline_dist, L) * 1e6 * (params$secpera)^params$m

plot(flowline_dist/1000, fric_sim / (1e6 * (params$secpera)^params$m), type = "l", 
    xlab = "Distance along flowline (km)", ylab = "Friction coefficient (Pa s^(1/3) m^(-1/3))")

## Surface elevation
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

## Initial thickness
se_grounded <- na.omit(surf_elev_mat[, 1]) # Use surface elevation in the year 2000 to initialise ice thickness
H_ini <- se_grounded - bed_sim[1:length(se_grounded)]
tail(H_ini)

## No thickness data on the ice shelf, so I'll just add on a shelf with gradually decreasing thickness
missing <- which(is.na(surf_elev_mat[, 1]))
thickness_at_gl <- - bed_sim[missing][1] * params$rho_w / params$rho_i
thickness_at_tail <- seq(from = thickness_at_gl, by = -1, length.out = length(missing)-1) #rep(400, length(missing)-1) #- bed_sim[missing] * rho_w / rho # worked out based on grounding line conditions

H_ini_all <- c(H_ini, thickness_at_gl, thickness_at_tail)

## Velocity
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_curr <- rowMeans(vel_mat[, 15:ncol(vel_mat)])

## Smooth the velocity out with loess
vel_curr_smooth <- loess(vel_curr ~ flowline_dist, span = 0.1)$fitted

## SMB data
smb_data_racmo <- qread(file = paste0(data_dir, "/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)

# ## Relaxation rate
# ## Change in surface elevation (used for relaxation)
# se_change <- rowMeans(surf_elev_mat[, 2:21] - surf_elev_mat[, 1:20]) # 2001 to 2020

# ## Change in shelf height
# shelf_height_change <- qread(file = paste0(data_dir, "/SMB/flowline_shelf_height_change.qs"))
# avg_shelf_height_change <- colMeans(shelf_height_change, na.rm = T)

# ## Fill in missing surface elevation change with shelf height change    
# missing <- which(is.na(se_change))
# se_change[missing] <- avg_shelf_height_change[missing] #-0.1 #tail(se_change[nonmissing], 1)
# se_change[is.na(se_change)] <- mean(avg_shelf_height_change, na.rm = T) ## fill in the rest of the missing values with the mean change in shelf height

# ## Now smooth this elevation change out with loess()
# se_change_smooth <- loess(se_change ~ flowline_dist, span = 0.08)$fitted

## Basal melt data (only available on shelves???)
if (use_basal_melt_data) {
  melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
  # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
  avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
  melt_nonmissing <- which(!is.na(avg_melt_rate))
  avg_melt_rate[1:(melt_nonmissing[1]-1)] <- -1 #seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
  avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) #mean(avg_melt_rate[melt_nonmissing])
  avg_melt_rate <- - avg_melt_rate # inverting this as eventually smb is calculated as smb - melt

} else {
  avg_melt_rate <- rep(1, J)
}

if (rerun_steady_state) {
    params$as <- smb_avg # surface accumulation rate (m/s)
    params$ab <- 0 # melt rate (m/s) -- no melt for now, just let the ice sheet grow
    steady_state <- solve_ssa_nl(domain = flowline_dist, 
                            bedrock = bed_sim, 
                            friction_coef = fric_sim, 
                            phys_params = params,
                            tol = 1e-02, #m/yr 
                            # years = 200,
                            steps_per_yr = 52, 
                            add_process_noise = F,
                            # thickness_bc = 3500,
                            ini_thickness = H_ini_all,
                            ini_velocity = vel_curr_smooth,
                            increase_hardness = T # this will increase ice rigidity parameter
                            # relax_rate = se_change_smooth,
                            # smb = smb_avg,
                            # basal_melt = 0 #basal_melt
                        )

    # qsave(steady_state, file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
   
} else {
    steady_state <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
}

## Then start ``thinning'' the ice by increasing melt rate to average present-day melt rate

### For this part I also add process noise
### Process noise parameters (for ice thickness)
exp_cov <- function(d, l) {
  return(exp(-3 * d / l))
}

ones <- rep(1, length(flowline_dist))
D <- rdist(flowline_dist)
l <- 50e3
R <- exp_cov(D, l)


# R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
L <- t(chol(R))
L <- as(L, "dgCMatrix")
process_noise_info <- list(corrmat_chol = L, length_scale = l)

relax_years <- 20 #00
relaxation <- solve_ssa_nl(domain = flowline_dist, 
                            bedrock = bed_sim, 
                            friction_coef = fric_sim, 
                            phys_params = params,
                            # tol = 1e-03, 
                            years = relax_years, #500,
                            steps_per_yr = 52, 
                            add_process_noise = T,
                            process_noise_info = process_noise_info,
                            ini_thickness = steady_state$current_thickness,
                            ini_velocity = steady_state$current_velocity,
                            # relax_rate = se_change_smooth,
                        )

# qsave(relaxation, file = paste0(data_dir, "training_data/steady_state/steady_state_relax_", data_date, ".qs"))

png(file = paste0("./plots/steady_state/relaxation_", data_date, ".png"), width = 800, height = 600)

par(mfrow = c(2,1))
matplot(flowline_dist/1000, steady_state$all_top_surface, type = "l", col = "grey", ylim = c(0, 1500),
    xlab = "Distance along flowline (km)", ylab = "Elevation (m)", main = paste0("Relaxation year ", relax_years))
lines(flowline_dist/1000, surf_elev_mat[, 21], col = "red", lty = 1)
# lines(flowline_dist/1000, steady_state$current_top_surface, col = "salmon")
# lines(flowline_dist/1000, relaxation$current_top_surface, col = "red")
legend("topright", legend = c("20 year post-steady state", "Present-day"), 
    col = c("grey", "red"), lty = 1, bty = "n")

matplot(flowline_dist/1000, relaxation$all_velocities, type = "l", col = "grey", ylim = c(0, 4000),
    xlab = "Distance along flowline (km)", ylab = "Velocity (m/a)")
lines(flowline_dist/1000, vel_curr_smooth, col = "red", lty = 1)
# lines(flowline_dist/1000, steady_state$current_velocity, col = "salmon")
# lines(flowline_dist/1000, relaxation$current_velocity, col = "red")
legend("topleft", legend = c("20 year post-steady state", "Present-day"), 
    col = c("grey", "red"), lty = 1, bty = "n")
dev.off()




