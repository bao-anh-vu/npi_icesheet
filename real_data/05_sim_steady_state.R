setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(dplyr)
library(fields)
library(Matrix)
library(matrixStats)
library(R.utils)
library(mvtnorm)
library(ggplot2)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library("sp")
# library(sf)

source("./source/create_params.R")
# source("./source/simulate_friction.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
source("./source/simulate_friction.R")
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
params$B <- 1.4 * 1e6 * params$secpera^params$m
params$A <- params$B^(-params$n)

## Flowline data
# flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## Simulate bed
# bed_sim <- create_bed(x = flowline_dist)
bed_sim <- qread(file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))

print("Simulating friction coefficient...")

## Simulate friction coefficient
# fric.sill <- 8e-5
# fric.nugget <- 0
# fric.range <- 10e3

# fric_sim <- simulate_friction2(
#     nsim = 1, domain = flowline_dist,
#     sill = fric.sill, nugget = fric.nugget,
#     range = fric.range
# ) 
L <- flowline_dist[J] - flowline_dist[1]
fric_sim <- create_fric_coef(flowline_dist, L) * 1e6 * (params$secpera)^params$m
# fric_val <- seq(0.01, 0.02, length.out = J) #* 1e6 * (params$secpera)^params$m
# fric_sim <- fric_val * 1e6 * (params$secpera)^params$m

png(file = paste0("./plots/steady_state/friction_coef_", data_date, ".png"), width = 800, height = 600)
plot(flowline_dist/1000, fric_sim / (1e6 * (params$secpera)^params$m), type = "l", 
    xlab = "Distance along flowline (km)", ylab = "Friction coefficient (Pa s^(1/3) m^(-1/3))")
dev.off()

## Grounding line position
gl_pos <- qread(file = paste0(data_dir, "grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind

## Surface elevation
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

## Initial thickness
se_grounded <- na.omit(rowMeans(surf_elev_mat)) # Use surface elevation in the year 2000 to initialise ice thickness
H_ini <- se_grounded - bed_sim[1:length(se_grounded)]
length_shelf <- J - length(se_grounded)

# se_gl <- se_grounded[gl_ind]
# se_shelf <- seq(from = se_gl, to = 100, length.out = length_shelf) # extend the same surface elevation at the GL to the ice shelf
# H_shelf <- - se_shelf * (params$rho_w / (params$rho_i - params$rho_w)) # thickness at grounding line based on flotation condition
# H_shelf <- se_shelf / (1 - params$rho_i / params$rho_w) # thickness at grounding line based on flotation condition
H_shelf <- - bed_sim[(gl_ind+1):J] * params$rho_w / params$rho_i - 100 # minus an offset to satisfy flotation condition 


# thickness_at_gl <- - bed_sim[gl_ind][1] * params$rho_w / params$rho_i
# H_shelf <- seq(thickness_at_gl - 1, 500, length.out = length_shelf)

H_ini_all <- c(H_ini, H_shelf)

z <- get_surface_elev(H = H_ini_all, b = bed_sim, z0 = 0, rho = params$rho_i, rho_w = params$rho_w, include_GL = TRUE)

png(file = paste0("./plots/steady_state/initial_thickness_", data_date, ".png"), width = 800, height = 600)
plot(flowline_dist/1000, z, type = "l")
abline(v = flowline_dist[length(se_grounded)]/1000, lty = 2, col = "red")
dev.off()

## Velocity
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
years <- ncol(vel_mat)
vel_curr <- rowMeans(vel_mat[, (years-5):years]) # average over the last 10 years

## Smooth the velocity out with gam()
vel_df <- data.frame(flowline_dist = flowline_dist, vel = vel_curr)
vel_data <- vel_df %>% filter(!is.na(vel))
vel_missing <- vel_df %>% filter(is.na(vel))
library(mgcv)
gam_fit <- gam(vel ~ s(flowline_dist, k = 20), data = vel_data)

# Predict at missing locations (works even if outside observed range)
vel_pred <- predict(
  gam_fit,
  newdata = data.frame(flowline_dist = vel_missing$flowline_dist)
)
vel_missing$vel <- vel_pred

vel_curr_smooth <- rbind(vel_data, vel_missing) %>%
  arrange(flowline_dist) %>%
  pull(vel)

## SMB data
smb_data_racmo <- qread(file = paste0(data_dir, "/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)
params$as <- smb_avg # surface accumulation rate (m/s)
params$ab <- 0 # melt rate (m/s) -- no melt for now, just let the ice sheet grow

if (rerun_steady_state) {
    print('Running model to steady state...')
    
    steady_state <- solve_ssa_nl(domain = flowline_dist, 
                            bedrock = bed_sim, 
                            friction_coef = fric_sim, 
                            phys_params = params,
                            tol = 1e-02, #m/yr 
                            # years = 100,
                            steps_per_yr = 100, 
                            add_process_noise = F,
                            # thickness_bc = 3500,
                            ini_thickness = H_ini_all,
                            ini_velocity = vel_curr_smooth
                        )

    qsave(steady_state, file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
   
} else {
    steady_state <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
}

png(file = paste0("./plots/steady_state/steady_state_", data_date, ".png"), width = 800, height = 600)

par(mfrow = c(2,1))
matplot(flowline_dist/1000, steady_state$current_top_surface, type = "l", col = "grey", ylim = c(0, 1500),
    xlab = "Distance along flowline (km)", ylab = "Elevation (m)", main = paste0("Relaxation year ", relax_years))
matlines(flowline_dist/1000, surf_elev_mat, col = "red", lty = 1)
# lines(flowline_dist/1000, relaxation$current_top_surface, col = "red")
legend("topright", legend = c(paste0("Simulated (", relax_years, "years post-steady)"), "Observed"), 
    col = c("grey", "red"), lty = 1, bty = "n")

matplot(flowline_dist/1000, steady_state$current_velocity, type = "l", col = "grey", ylim = c(0, 4000),
    xlab = "Distance along flowline (km)", ylab = "Velocity (m/a)")
matlines(flowline_dist/1000, vel_mat, col = "red", lty = 1)
# lines(flowline_dist/1000, steady_state$current_velocity, col = "salmon")
# lines(flowline_dist/1000, relaxation$current_velocity, col = "red")
legend("topleft", legend = c("Simulated (20 year post-steady)", "Observed"), 
    col = c("grey", "red"), lty = 1, bty = "n")
dev.off()

# ## Then start ``thinning'' the ice by using present-day average melt rate
# ## Basal melt data (only available on shelves???)
# if (use_basal_melt_data) {
#   melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
#   # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
#   avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
#   melt_nonmissing <- which(!is.na(avg_melt_rate))
#   avg_melt_rate[1:(melt_nonmissing[1]-1)] <- -1 #seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
#   avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) #mean(avg_melt_rate[melt_nonmissing])
#   avg_melt_rate <- - avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
# } else {
#   avg_melt_rate <- rep(0, J)
# }
# params$ab <- avg_melt_rate # melt rate (m/s)

# ## Plot average melt rate
# png(file = paste0("./plots/steady_state/avg_melt_rate_", data_date, ".png"), width = 800, height = 600)
# plot(flowline_dist/1000, params$ab, type = "l", 
#     xlab = "Distance along flowline (km)", ylab = "Basal melt rate (m/a)")
# abline(v = flowline_dist[length(se_grounded)]/1000, lty = 2, col = "red")
# dev.off()

# ### Also reduce ice rigidity
# params$B <- 1.2 * 1e6 * params$secpera^params$m

# ### For this part I also add process noise
# ### Process noise parameters (for ice thickness)
# exp_cov <- function(d, l) {
#   return(exp(-3 * d / l))
# }

# ones <- rep(1, length(flowline_dist))
# D <- rdist(flowline_dist)
# l <- 50e3
# R <- exp_cov(D, l)

# # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
# L <- t(chol(R))
# L <- as(L, "dgCMatrix")
# process_noise_info <- list(corrmat_chol = L, length_scale = l)

# relax_years <- 10
# cat("Running model for ", relax_years, "years post-steady state...")

# relaxation <- solve_ssa_nl(domain = flowline_dist, 
#                             bedrock = bed_sim, 
#                             friction_coef = fric_sim, 
#                             phys_params = params,
#                             # tol = 1e-03, 
#                             years = relax_years, #500,
#                             steps_per_yr = 100, 
#                             add_process_noise = F,
#                             process_noise_info = process_noise_info,
#                             ini_thickness = steady_state$current_thickness,
#                             ini_velocity = steady_state$current_velocity#,
#                             # use_relaxation = T,
#                             # observed_thickness = H_ini_all
#                             # relax_rate = relax_rate
#                         )

# qsave(relaxation, file = paste0(data_dir, "training_data/steady_state/steady_state_relax_", data_date, ".qs"))

# png(file = paste0("./plots/steady_state/relaxation_", data_date, ".png"), width = 800, height = 600)

# par(mfrow = c(2,1))
# matplot(flowline_dist/1000, relaxation$all_top_surface, type = "l", col = "grey", ylim = c(0, 1500),
#     xlab = "Distance along flowline (km)", ylab = "Elevation (m)", main = paste0("Relaxation year ", relax_years))
# matlines(flowline_dist/1000, surf_elev_mat, col = "red", lty = 1)
# lines(flowline_dist/1000, steady_state$current_top_surface, col = "salmon")
# # lines(flowline_dist/1000, relaxation$current_top_surface, col = "red")
# legend("topright", legend = c(paste0("Simulated (", relax_years, "years post-steady)"), "Observed"), 
#     col = c("grey", "red"), lty = 1, bty = "n")

# matplot(flowline_dist/1000, relaxation$all_velocities, type = "l", col = "grey", ylim = c(0, 4000),
#     xlab = "Distance along flowline (km)", ylab = "Velocity (m/a)")
# matlines(flowline_dist/1000, vel_mat, col = "red", lty = 1)
# lines(flowline_dist/1000, steady_state$current_velocity, col = "salmon")
# # lines(flowline_dist/1000, relaxation$current_velocity, col = "red")
# legend("topleft", legend = c("Simulated (20 year post-steady)", "Observed"), 
#     col = c("grey", "red"), lty = 1, bty = "n")
# dev.off()




