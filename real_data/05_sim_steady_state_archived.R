# Simulate steady state
setwd("~/SSA_model/CNN/real_data/")

library(dplyr)
library(ggplot2)
library(RColorBrewer)
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
library(FRK)
library(qs)
# library(sf)

source("./source/create_params.R")
source("./source/cond_sim_gp.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/sim_steady_state.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
source("./source/fit_basis.R")

data_dir <- "./data/"
data_date <- "20241111" #"20241103"

## Flags
simulate_steady_state <- T
# smooth_bed <- F
save_steady_state <- F

## Flowline data
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
J <- nrow(flowline) # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## 1. Read bed data
bed_sim <- qread(file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))

## Simulate friction coefficient
print("Simulating friction coefficient...")

secpera <- 31556926 #seconds per annum
n <- 3.0 # exponent in Glen's flow law
m <- 1/n # friction exponent
x <- flowline_dist
L <- flowline_dist[J] - flowline_dist[1]
fric_sim <- create_fric_coef(x, L) * 1e6 * (secpera)^m

## Grounding line pos
# gl_pos <- readRDS(file = "./data/grounding_line/gl_pos.rds")

## Surface elevation
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

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

plot(flowline_dist/1000, se_change, type = "l", xlab = "Domain (km)", ylab = "Elevation change (m/a)")
lines(flowline_dist/1000, se_change_smooth, col = "cyan")
abline(v = flowline_dist[missing[1]]/1000, col = "red", lty = 2) # grounding line

## Initial thickness
se_grounded <- na.omit(surf_elev_mat[, 21]) # Use surface elevation at the final time to initialise ice thickness
H_ini <- se_grounded - bed_sim[1:length(se_grounded)]
tail(H_ini)
missing <- which(is.na(surf_elev_mat[, 21]))
rho <- 910.0
rho_w <- 1028.0
thickness_at_gl <- - bed_sim[missing][1] * rho_w / rho
thickness_at_tail <- rep(400, length(missing)-1) #- bed_sim[missing] * rho_w / rho # worked out based on grounding line conditions

H_ini_new <- c(H_ini, thickness_at_gl, thickness_at_tail) #+ offset
# H_ini_new <- H_ini_new + offset

## Velocity
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_curr <- rowMeans(vel_mat[, 15:ncol(vel_mat)])

## Smooth the velocity out with loess
vel_curr_smooth <- loess(vel_curr ~ flowline_dist, span = 0.1)$fitted

# png(paste0("./plots/steady_state/vel_curr.png"), width = 800, height = 600)
# plot(vel_curr, type = "l")
# lines(fitted(vel_curr_smooth), col = "red")
# dev.off()

# png(paste0("./plots/steady_state/H_ini.png"), width = 1000, height = 500)
# plot(H_ini_new, type = "l")
# dev.off()

## SMB data
smb_data_racmo <- qread(file = paste0(data_dir, "/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_shelf_data <- qread(file = paste0(data_dir, "/SMB/flowline_shelf_smb.qs")) ## from 1992 to 2017
melt_data <- qread(file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))

smb_avg <- colMeans(smb_data_racmo, na.rm = T)
smb_shelf_avg <- colMeans(smb_shelf_data, na.rm = T)
melt_avg <- colMeans(melt_data, na.rm = T)
melt_avg[is.na(melt_avg)] <- 0

## Run model to steady state based on bed_sim
if (simulate_steady_state) {
        print("Simulating steady state...")
        # steady_state <- sim_steady_state(domain = flowline_dist,
        #                         years = 100,
        #                         bedrock = bed_sim,
        #                         friction = fric_sim,
        #                         # ini_thickness = H_ini_new,
        #                         # ini_velocity = vel_curr_smooth$fitted,
        #                         # ini_velocity = 0.001 / secpera * flowline_dist,
        #                         # relax_rate = se_change_smooth,
        #                         smb = smb_avg,
        #                         basal_melt = 0 #melt_avg
        #                         )

        steady_state <- solve_ssa_nl(domain = flowline_dist, 
                                bedrock = bed_sim, 
                                friction_coef = fric_sim, 
                                tol = 1e-03, 
                                #     years = 50,
                                steps_per_yr = 100, 
                                #     save_model_output = TRUE, 
                                perturb_hardness = FALSE, 
                                add_process_noise = F,
                                thickness_bc = 3000,
                                # ini_thickness = H_ini_new,
                                ini_velocity = vel_curr_smooth,
                                # relax_rate = se_change_smooth,
                                smb = smb_avg,
                                basal_melt = 0 #basal_melt
                            )

        if (save_steady_state) {
                # saveRDS(steady_state, file = paste0("./data/training_data/steady_state_", date, ".rds"))
                qsave(steady_state, file = paste0(data_dir, "/steady_state_", data_date, ".qs"))

        }

        # lines(flowline_dist, bedrock)
} else {
        steady_state <- qread(file = paste0(data_dir, "/steady_state_", data_date, ".qs"))
}

## Plot the ice sheet at steady state

domain <- steady_state$domain
thickness <- steady_state$current_thickness
bedrock <- steady_state$bedrock
top_surface_elev <- steady_state$top_surface
bottom_surface_elev <- steady_state$bottom_surface

png(paste0("./plots/steady_state_", data_date, ".png"), width = 800, height = 800)
plot(domain, top_surface_elev, type = 'l', ylim = c(-2000, 2000))
lines(domain, bottom_surface_elev)
lines(domain, bedrock)
lines(domain, surf_elev_mat[, 21])
dev.off()

gl <- steady_state$grounding_line
