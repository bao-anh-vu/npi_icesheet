## Test to see if a steady state can be simulated based on current SMB 
## but with the bed in the simulation study

setwd("~/SSA_model/CNN/real_data/")

library(dplyr)
library(ggplot2)
# library(RColorBrewer)
library(fields)
# library(parallel)
library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
library(R.utils)
# library("sp")
# library("tidyr")
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
# library(FRK)
library(qs)
# library(sf)

source("./source/create_params.R")
# source("./source/cond_sim_gp.R")
# source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
# source("./source/sim_steady_state.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
# source("./source/fit_basis.R")

data_dir <- "./data/"
data_date <- "20241111" #"20241103"

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
secpera <- 31556926 #seconds per annum
n <- 3.0 # exponent in Glen's flow law
m <- 1/n # friction exponent
L <- flowline_dist[J] - flowline_dist[1]
fric_sim <- create_fric_coef(flowline_dist, L) * 1e6 * (secpera)^m

plot(flowline_dist/1000, fric_sim / (1e6 * (secpera)^m), type = "l", 
    xlab = "Distance along flowline (km)", ylab = "Friction coefficient (Pa s^(1/3) m^(-1/3))")

## Surface elevation
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

## Initial thickness
se_grounded <- na.omit(surf_elev_mat[, 1]) # Use surface elevation in the year 2000 to initialise ice thickness
H_ini <- se_grounded - bed_sim[1:length(se_grounded)]
tail(H_ini)
missing <- which(is.na(surf_elev_mat[, 21]))
rho <- 910.0
rho_w <- 1028.0
thickness_at_gl <- - bed_sim[missing][1] * rho_w / rho
thickness_at_tail <- seq(from = thickness_at_gl, by = -1, length.out = length(missing)-1) #rep(400, length(missing)-1) #- bed_sim[missing] * rho_w / rho # worked out based on grounding line conditions

H_ini_new <- c(H_ini, thickness_at_gl, thickness_at_tail) #+ offset

## Velocity
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_curr <- rowMeans(vel_mat[, 15:ncol(vel_mat)])

## Smooth the velocity out with loess
vel_curr_smooth <- loess(vel_curr ~ flowline_dist, span = 0.1)$fitted

## SMB data
smb_data_racmo <- qread(file = paste0(data_dir, "/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)

## Relaxation rate
## Change in surface elevation (used for relaxation)
se_change <- rowMeans(surf_elev_mat[, 2:21] - surf_elev_mat[, 1:20]) # 2001 to 2020

## Change in shelf height
shelf_height_change <- qread(file = paste0(data_dir, "/SMB/flowline_shelf_height_change.qs"))
avg_shelf_height_change <- colMeans(shelf_height_change, na.rm = T)

## Fill in missing surface elevation change with shelf height change    
missing <- which(is.na(se_change))
se_change[missing] <- avg_shelf_height_change[missing] #-0.1 #tail(se_change[nonmissing], 1)
se_change[is.na(se_change)] <- mean(avg_shelf_height_change, na.rm = T) ## fill in the rest of the missing values with the mean change in shelf height

## Now smooth this elevation change out with loess()
se_change_smooth <- loess(se_change ~ flowline_dist, span = 0.08)$fitted

steady_state <- solve_ssa_nl(domain = flowline_dist, 
                            bedrock = bed_sim, 
                            friction_coef = fric_sim, 
                            tol = 1e-03, 
                            # years = 200,
                            steps_per_yr = 100, 
                            add_process_noise = F,
                            # thickness_bc = 3500,
                            ini_thickness = H_ini_new,
                            ini_velocity = vel_curr_smooth,
                            perturb_hardness = T, # this will increase ice rigidity parameter
                            # relax_rate = se_change_smooth,
                            smb = smb_avg,
                            basal_melt = 0 #basal_melt
                        )

# qsave(steady_state, file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
