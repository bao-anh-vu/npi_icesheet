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
library(mgcv)

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
smooth_bed <- F
set.seed(2025)
## Flowline data
# flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
J <- nrow(flowline)
# J <- 2001 # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## 1. Read bed data
# bed_old <- readRDS(file = "./data/bed_elev_nearest.rds")
# bed <- readRDS(file = "./data/bedmap_obs.rds")
# bed_sd <- unlist(readRDS(file = "./data/bedmap_sd.rds"))
bed <- qread(file = paste0(data_dir, "bedmap_obs.qs"))
bed_sd <- unlist(qread(file = paste0(data_dir, "bedmap_sd.qs")))
bed <- bed[1:J] 
bed_sd <- bed_sd[1:J]

bed_df <- data.frame(ind = 1:J, loc = flowline_dist, bed_elev = bed, bed_sd = bed_sd)

# ## Calculate bed elevation at GL based on flotation condition
gl_pos <- qread(file = paste0(data_dir, "grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind

# ## Condition at GL: bed = surface elevation / (rho_w/rho_i + 1)
# params <- list(
#   rho_i = 917.0, # ice density
#   rho_w = 1028.0 # sea water density
# )
# surf_elev_mat <- qread(paste0(data_dir, "surface_elev/surf_elev_mat.qs"))
# se_gl <- mean(surf_elev_mat[gl_ind, ]) # use the time-averaged surface observation at GL

# bed_gl <- se_gl / (1 - params$rho_w / params$rho_i)
# bed_obs_df$bed_elev[gl_ind] <- bed_gl
# bed_obs_df$bed_sd[gl_ind] <- 0 

# Construct bed based on observations
bed_obs_df <- bed_df %>% filter(!is.na(bed_elev) & !is.na(bed_sd)) #%>% pull(ind)

# Choose some bed observations; leave the rest for validation
# Divide domain into intervals of 10km
interval_length <- 5e3 # m
# # num_intervals <- floor((max(flowline_dist) - min(flowline_dist)) / interval_length)
cut_points <- seq(from = min(flowline_dist), to = max(flowline_dist), by = interval_length)
# intervals <- cut(flowline_dist, breaks = cut_points, include.lowest = TRUE, right = FALSE)

# Sort bed observations into intervals
bed_obs_df$interval <- cut(bed_obs_df$loc, breaks = cut_points, include.lowest = TRUE, right = FALSE)
intervals <- unique(bed_obs_df$interval)

# Alternate between 0 and 1 for the length of the 'intervals' vector
chosen_int <- rep(c(0, 1), length.out = length(intervals))

# Select every other interval
intervals_df <- data.frame(interval = intervals, chosen = chosen_int)
# selected_bed_obs <- bed_obs_df %>% filter(interval %in% selected_intervals

# Now label whether a bed observation is chosen based on its interval
bed_obs_df <- bed_obs_df %>% left_join(intervals_df, by = "interval")

## set.seed(2025)
## chosen_bed_ind <- sample(common_avail_ind, 50, replace = F)

## if index is in chosen_bed_ind, then chosen = 1, else 0
## bed_obs_df$chosen <- ifelse(bed_obs_df$ind %in% chosen_bed_ind, 1, 0)
### bed_obs_df$chosen[gl_ind] <- 1 # also choose the bed "observation" at GL

qsave(bed_obs_df, file = paste0(data_dir, "/bed_obs_df.qs"))

bed_obs_chosen <- bed_obs_df %>% filter(chosen == 1) %>% na.omit()


# bed.sill <- 5e3
# bed.range <- 50e3
# bed.nugget <- 0 #200
n_sims <- 10
bed_sim_output <- simulate_bed(n_sims, domain = flowline_dist, 
                        obs_locations = bed_obs_chosen$ind, 
                        obs = bed_obs_chosen$bed_elev) #, 
                        # sill = bed.sill, nugget = bed.nugget,
                        # range = bed.range)
# bed_sim_output <- simulate_bed(1, domain = flowline_dist, 
#                             obs_locations = bed_obs_chosen$ind, 
#                             obs = bed_obs_chosen$bed_elev) 

## Turn bed_sim into a data frame for plotting
bed_sim_df <- data.frame(x = flowline_dist, y = bed_sim_output$sims)
colnames(bed_sim_df) <- c("x", paste0("sim", 1:n_sims))
bed_mean <- bed_sim_output$mean

## Plot the bed simulations
bed_sim_plot <- bed_sim_df %>% ggplot(aes(x = loc, y = bed_elev)) +
        geom_line(data = bed_sim_df, aes(x = x, y = sim1), col = "grey40") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim2), col = "grey60") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim3), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim4), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim5), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim6), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim7), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim8), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim9), col = "grey80") +
        geom_line(data = bed_sim_df, aes(x = x, y = sim10), col = "grey80") +
        geom_line(data = data.frame(x = flowline_dist, y = bed_mean), aes(x = x, y = y), col = "red", lwd = 1) +
        geom_point(data = bed_obs_df, aes(col = factor(chosen)), size = 3) + #, shape = 21, size = 3) +
        geom_vline(xintercept = flowline_dist[gl_ind], lty = 2) +
        xlim(c(0, 200e3)) + 
        theme_bw() +
        theme(text = element_text(size = 20))

png(paste0("./plots/bed/bed_simulations.png"), width = 1000, height = 500)
print(bed_sim_plot)
dev.off()

## Choose one bed simulation for the steady state
bed_sim <- bed_sim_output$sims[, 1] 

## Smooth tbe bed simulation with loess()
if (smooth_bed) {
        bed_loess <- loess(bed ~ x, data = data.frame(bed = bed_sim, x = flowline_dist), span = 0.01)
        bed_sim <- fitted(bed_loess)
}

## Alternative: use GP regression
## Add "artificial" observations at the first and last points to constrain the bed
## Otherwise the bed will blow up at the boundaries
# bed_obs_append <- c(bed_obs_chosen$bed_elev[1], bed_obs_chosen$bed_elev, 
#                         bed_obs_chosen$bed_elev[length(bed_obs_chosen$bed_elev)])
# bed_obs_loc_append <- c(flowline_dist[1], bed_obs_chosen$loc, 
#                         flowline_dist[length(flowline_dist)])
# bed_pred_loc <- bed_obs_df$loc[bed_obs_df$chosen == 0]

# out <- cond_sim_gp(10, x_test = bed_pred_loc, 
#                         x_train = bed_obs_chosen$loc,
#                         obs = bed_obs_chosen$bed_elev, 
#                         obs_sd = rep(0, length(bed_obs_chosen$chosen == 0))) # use zero sd

# out2 <- cond_sim_gp(10, x_test = bed_pred_loc, 
#                         x_train = bed_obs_chosen$loc,
#                         obs = bed_obs_chosen$bed_elev, 
#                         obs_sd = bed_obs_chosen$bed_sd)

# matplot(out$sims$x, out$sims[, -(1:2)], type = "l", col = "grey")
# points(bed_obs_chosen$loc[bed_obs_chosen$chosen == 1], bed_obs_chosen$bed_elev[bed_obs_chosen$chosen == 1])

# matplot(out$sims$x, out2$sims[, -(1:2)], type = "l", col = "grey")
# points(bed_obs_chosen$loc, bed_obs_chosen$bed_elev)

# bed_sim2 <- out$sims[, 2]
# bed_sim3 <- out$sims[, 3]
# bed_sim4 <- out$sims[, 4]
bed_sim_plot <- bed_obs_df %>% ggplot(aes(x = loc, y = bed_elev)) +
        # geom_line() + 
        geom_point(aes(col = factor(chosen)), size = 3) + #, shape = 21, size = 3) +
        geom_vline(xintercept = flowline_dist[gl_ind], lty = 2) +
        geom_line(data = data.frame(x = flowline_dist, y = bed_sim), aes(x = x, y = y), col = "grey60") +
        xlim(c(0, 200e3)) + 
        theme_bw() +
        theme(text = element_text(size = 20)) 

png(paste0("./plots/bed/bed_sim.png"), width = 1000, height = 500)
print(bed_sim_plot)
dev.off()

bed_obs_chosen <- bed_obs_df %>% filter(chosen == 1) 

png(paste0("./plots/bed/bed_obs_chosen.png"), width = 1000, height = 500)
bed_obs_chosen %>% ggplot(aes(x = loc, y = bed_elev)) +
        geom_point(alpha = 0.5, size = 3) + #, shape = 21, size = 3) +
        geom_vline(xintercept = flowline_dist[gl_ind], lty = 2) +
        geom_line(data = data.frame(x = flowline_dist, y = bed_sim), aes(x = x, y = y), col = "grey60") +
        # geom_line(data = data.frame(x = flowline_dist, y = fitted(bed_loess)), aes(x = x, y = y), col = "red") +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim2), aes(x = x, y = y), col = "blue") +
        # # xlim(2e+05, 3e+05) +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim3), aes(x = x, y = y), col = "green") +
        # geom_line(data = data.frame(x = flowline_dist, y = bed_sim4), aes(x = x, y = y), col = "purple") +
        theme_bw() +
        theme(text = element_text(size = 20)) 
dev.off()

## Save
qsave(bed_sim, file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))
