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
library(fields)
# library("tidyr")
library(dplyr)
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(FRK)
library(qs)
# library(sf)
source("./source/create_params.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/sim_steady_state.R")
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
source("./source/fit_basis.R")
source("./source/create_params.R")

data_dir <- "./data/"
data_date <- "20241103"

## Flags
simulate_steady_state <- T
use_basis_funs <- F
save_steady_state <- T

## Flowline data
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
J <- 2001 # number of grid points
flowline <- flowline[1:J, ]
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

bed_sim <- simulate_bed(1, domain = flowline_dist, 
                        obs_locations = chosen_bed_df$ind, 
                        obs = chosen_bed_df$bed_elev, obs_sd = chosen_bed_df$bed_sd)


png(paste0("./plots/bed/bed_sim.png"), width = 1000, height = 500)
bed_obs_df %>% ggplot(aes(x = loc, y = bed_elev)) +
        geom_line() + 
        geom_point(aes(col = factor(chosen))) +
        geom_vline(xintercept = flowline_dist[gl_ind], lty = 2) +
        geom_line(data = data.frame(x = flowline_dist, y = bed_sim), aes(x = x, y = y), col = "gray") +
        theme_bw()
dev.off()


print("Simulating friction coefficient...")

fric.sill <- 8e-5
fric.nugget <- 0
fric.range <- 10e3

fric_sim <- simulate_friction2(
nsim = 1, domain = flowline_dist,
sill = fric.sill, nugget = fric.nugget,
range = fric.range
) 

secpera <- 31556926 #seconds per annum
n <- 3.0 # exponent in Glen's flow law
m <- 1/n # friction exponent
x <- flowline_dist
L <- flowline_dist[J] - flowline_dist[1]
fric_sim <- create_fric_coef(x, L) * 1e6 * (secpera)^m

# png(paste0("./plots/friction/fric_steady_state.png"), width = 800, height = 800)
# plot(flowline_dist, fric_sim, type = "l")
# dev.off()

# H0 <- 1000
# L <- flowline_dist[J] - flowline_dist[1]
# x <- flowline_dist
# H_ini <- H0 - (H0 - 0)/(L - 0) * x
# # elev_at_tail <- 100
# # H_ini <- H0 - (H0 - elev_at_tail)/(L^2) * x^2

# plot(bed_sim, type = "l", ylim = c(-2000, 2000))
# lines(H_ini)


## Run model to steady state based on bed_sim
if (simulate_steady_state) {
        print("Simulating steady state...")
        steady_state <- sim_steady_state(domain = flowline_dist,
                                bedrock = bed_sim,
                                friction = fric_sim
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

png(paste0("./plots/steady_state.png"), width = 800, height = 800)
plot(domain, top_surface_elev, type = 'l', ylim = c(-2000, 2000))
lines(domain, bottom_surface_elev)
lines(domain, bedrock)
dev.off()

gl <- steady_state$grounding_line
