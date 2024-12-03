## Simulate bed for steady state# Simulate steady state
setwd("~/SSA_model/CNN/real_data/")

library(dplyr)
library(ggplot2)
# library(RColorBrewer)
library(fields)
# library(parallel)
# library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library(R.utils)
# library("sp")
# library("tidyr")
# library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
# library(FRK)
library(qs)

# library(sf)
# source("./source/create_params.R")
# source("./source/cond_sim_gp.R")
source("./source/simulate_bed.R")
# source("./source/simulate_friction.R")
# source("./source/sim_steady_state.R")
# source("./source/solve_ssa_nl.R")
# source("./source/solve_velocity_azm.R")
# source("./source/solve_thickness.R")
source("./source/surface_elev.R")
# source("./source/fit_basis.R")

data_dir <- "./data/"
data_date <- "20241103"

## Flags
smooth_bed <- T

## Flowline data
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
J <- 2001 # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## 1. Read bed data
# bed_old <- readRDS(file = "./data/bed_elev_nearest.rds")
bed <- readRDS(file = "./data/bedmap_obs.rds")
bed_sd <- unlist(readRDS(file = "./data/bedmap_sd.rds"))
bed <- bed[1:J] 
bed_sd <- bed_sd[1:J]

bed_obs_df <- data.frame(ind = 1:J, loc = flowline_dist, bed_elev = bed, bed_sd = bed_sd)
bed_avail_ind <- which(!is.na(bed))
bed_sd_avail_ind <- which(!is.na(bed_sd))
common_avail_ind <- intersect(bed_avail_ind, bed_sd_avail_ind)

## Choose the first and last 25 observations, then randomly choose another 50 in between
chosen_bed_ind_head <- common_avail_ind[1:25] #Choose the first 50 bed observations
chosen_bed_ind_tail <- common_avail_ind[(length(common_avail_ind) - 24):length(common_avail_ind)] #Choose the first 50 bed observations
chosen_bed_ind <- c(chosen_bed_ind_head, chosen_bed_ind_tail)
set.seed(2024)
chosen_bed_ind_mid <- sample(setdiff(common_avail_ind, chosen_bed_ind), 50) #then randomly select another 50 bed observations after
chosen_bed_ind <- c(chosen_bed_ind_head, chosen_bed_ind_mid, chosen_bed_ind_tail)

# if index is in chosen_bed_ind, then chosen = 1, else 0
bed_obs_df$chosen <- ifelse(bed_obs_df$ind %in% chosen_bed_ind, 1, 0)

qsave(bed_obs_df, file = paste0(data_dir, "/bed_obs_df.qs"))

## Mark GL position
gl_pos <- readRDS(file = paste0(data_dir, "/grounding_line/gl_pos.rds"))
delta <- 500
pts_near_gl <- flowline %>% filter(
        x >= (gl_pos[1] - delta) & x <= (gl_pos[1] + delta),
        y >= (gl_pos[2] - delta) & y <= (gl_pos[2] + delta)
        ) %>%
        mutate(dist = sqrt((x - gl_pos[1])^2 + (y - gl_pos[2])^2)) %>%
        # arrange(dist) %>%
        slice_min(dist, n = 1) #%>%

gl_ind <- which(flowline$x == pts_near_gl$x)

chosen_bed_df <- bed_obs_df %>% filter(chosen == 1) %>% na.omit()

# bed_sim <- simulate_bed(1, domain = flowline_dist, 
#                         obs_locations = chosen_bed_df$ind, 
#                         obs = chosen_bed_df$bed_elev, 
#                         obs_sd = chosen_bed_df$bed_sd)

bed_sim_output <- simulate_bed(1, domain = flowline_dist, 
                            obs_locations = chosen_bed_df$ind, 
                            obs = chosen_bed_df$bed_elev) 
bed_sim <- bed_sim_output$sims
bed_mean <- bed_sim_output$mean

## Smooth tbe bed simulation with loess()
if (smooth_bed) {
        bed_loess <- loess(bed ~ x, data = data.frame(bed = bed_sim, x = flowline_dist), span = 0.01)
        bed_sim <- fitted(bed_loess)
}

## Alternative: use GP regression
## Add "artificial" observations at the first and last points to constrain the bed
## Otherwise the bed will blow up at the boundaries
# bed_obs_append <- c(chosen_bed_df$bed_elev[1], chosen_bed_df$bed_elev, 
#                         chosen_bed_df$bed_elev[length(chosen_bed_df$bed_elev)])
# bed_obs_loc_append <- c(flowline_dist[1], chosen_bed_df$loc, 
#                         flowline_dist[length(flowline_dist)])
# bed_pred_loc <- bed_obs_df$loc[bed_obs_df$chosen == 0]

# out <- cond_sim_gp(10, x_test = bed_pred_loc, 
#                         x_train = chosen_bed_df$loc,
#                         obs = chosen_bed_df$bed_elev, 
#                         obs_sd = rep(0, length(chosen_bed_df$chosen == 0))) # use zero sd

# out2 <- cond_sim_gp(10, x_test = bed_pred_loc, 
#                         x_train = chosen_bed_df$loc,
#                         obs = chosen_bed_df$bed_elev, 
#                         obs_sd = chosen_bed_df$bed_sd)

# matplot(out$sims$x, out$sims[, -(1:2)], type = "l", col = "grey")
# points(chosen_bed_df$loc[chosen_bed_df$chosen == 1], chosen_bed_df$bed_elev[chosen_bed_df$chosen == 1])

# matplot(out$sims$x, out2$sims[, -(1:2)], type = "l", col = "grey")
# points(chosen_bed_df$loc, chosen_bed_df$bed_elev)

# bed_sim2 <- out$sims[, 2]
# bed_sim3 <- out$sims[, 3]
# bed_sim4 <- out$sims[, 4]

png(paste0("./plots/bed/bed_sim.png"), width = 1000, height = 500)
bed_obs_df %>% ggplot(aes(x = loc, y = bed_elev)) +
        geom_line() + 
        geom_point(aes(col = factor(chosen))) +
        geom_vline(xintercept = flowline_dist[gl_ind], lty = 2) +
        geom_line(data = data.frame(x = flowline_dist, y = bed_sim), aes(x = x, y = y), col = "gray") +
        # geom_line(data = data.frame(x = flowline_dist, y = fitted(bed_loess)), aes(x = x, y = y), col = "red") +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim2), aes(x = x, y = y), col = "blue") +
        # # xlim(2e+05, 3e+05) +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim3), aes(x = x, y = y), col = "green") +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim4), aes(x = x, y = y), col = "purple") +
        theme_bw()
dev.off()

## Save
qsave(bed_sim, file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))
