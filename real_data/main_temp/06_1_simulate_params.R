# NOTE: need to either smooth out the bed simulations or
# find a way to modify the parameters in the conditional sim so that the bed isn't so noisy

## Generate 1000 simulations

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library(parallel)
# library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library(R.utils)
# library("sp")
library(fields)
# library("tidyr")
library(dplyr)
# library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(abind)
library(ggplot2)
library(gridExtra)
# library(FRK)
library(qs)
# library(mgcv)

# source("./source/sim_obs.R")
# source("./source/cond_sim_gp.R")
# source("./source/get_ini_thickness.R")
# source("./source/process_sim_results.R")
# source("./source/fit_basis.R")
# source("./source/surface_elev.R")
# source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
# source("./source/azm_cond_sim.R")

## Some flags
regenerate_sims <- T
# refit_basis <- T
save_sims <- T
# log_transform <- T
use_basal_melt_data <- T
constrain_gl <- F

## Directory for training data
train_data_dir <- "./data/training_data"

## Presets
data_date <- "20241111"
N <- 1000 # number of simulations per set
# set <- 1 #commandArgs(trailingOnly = TRUE)
sets <- 11:20 #50 #:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])
# warmup <- 0
# years <- 11 # number of years data is collected (not including initial condition)

## Physical params
# params <- list(
#   secpera = 31556926, # seconds per annum
#   n = 3.0, # exponent in Glen's flow law
#   rho_i = 917.0, # ice density
#   rho_w = 1028.0, # sea water density
#   g = 9.81 # gravity constant
#   # A = 4.227e-25, #1.4579e-25, # flow rate parameter
# )

# params$m <- 1 / params$n
# params$B <- 0.6 * 1e6 * params$secpera^params$m
# params$A <- params$B^(-params$n)

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# 0. Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# 0. Read bed observations
# bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
# bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]
bed_obs_df <- qread(file = paste0("./data/bedmap/bed_obs_df_all.qs"))
bed_obs_chosen <- bed_obs_df

# if (constrain_gl) {
#   ## Impose an exta condition at GL: bed = surface elevation / (rho_w/rho_i + 1)

#   steady_se <- ssa_steady$current_top_surface
#   steady_gl <- ssa_steady$grounding_line[length(ssa_steady$grounding_line)]
#   gl_ind <- sum(domain/1e3 < steady_gl)
#   bed_gl <- steady_se[gl_ind] / (1 - params$rho_w / params$rho_i)
#   bed_gl_df <- data.frame(ind = gl_ind, loc = domain[gl_ind], bed_elev = bed_gl, bed_sd = 0, interval = NA, chosen = 1)
#   bed_obs_df <- rbind(bed_obs_df, bed_gl_df) %>% arrange(ind)
# }
# bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]
# bed_obs_val <- bed_obs_df[bed_obs_df$chosen == 0, ]

## Bedmachine data to compare
bedmachine_data <- qread(paste0("./data/bedmachine/flowline_bedmachine.qs"))
bedmachine <- bedmachine_data$bed_avg

## Observed grounding line position
gl_pos <- qread(file = paste0("./data/grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind

# ## SMB data
# smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
# smb_avg <- colMeans(smb_data_racmo, na.rm = T)
# params$as <- smb_avg # surface accumulation rate (m/s)

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

t1 <- proc.time()
# if (regenerate_sims) {

flags <- c()

bed_sim_list <- list()
fric_sim_list <- list()  

## Bed prior from GP fit to BedMap data
bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))

for (i in 1:length(sets)) {
  set <- sets[i]
  setf <- formatC(sets[i], width = 2, flag = "0")
    
  cat("Generating set", set, "\n")

  ## Simulate bed and friction
  print("Simulating friction coefficient...")

  fric_sims <- simulate_friction2(
    nsim = N, domain = domain
  ) 

  ## Simulate beds
  print("Simulating beds...")

  # bed.sill <- 10e3
  # bed.range <- 50e3
  # bed.nugget <- 0 #200
  # bed_sim_output <- simulate_bed(N, domain = domain,
  #                         obs_locations = bed_obs_chosen$ind,
  #                         obs = bed_obs_chosen$bed_elev) #,

  # bed_sims <- bed_sim_output$sims
  # bed_mean <- bed_sim_output$mean


  L <- t(chol(bed_prior$cov))
  mean_mat <- matrix(rep(bed_prior$mean, N), nrow = nrow(bed_prior$mean), ncol = N)
  # L <- t(chol(bed_prior$cov))
  u_mat <- matrix(rnorm(nrow(L) * N), nrow = nrow(L), ncol = N)
  # mean_mat <- matrix(rep(bed_prior$mean, N), nrow = nrow(bed_prior$mean), ncol = N)
  bed_sims <- mean_mat + L %*% u_mat # rnorm(nrow(L) * N)

  # png(file = paste0("./plots/cnn/input/bed_sims_", setf, "_", data_date, ".png"), width = 800, height = 500)
  # matplot(domain, bed_sims[, 5:6], type = "l", col = "salmon", ylab = "Bedrock Elev (m)", xlab = "Distance along flowline (m)")
  # lines(domain, bedmachine, col = "blue", lwd = 2)
  # points(bed_obs_df$loc, bed_obs_df$bed_elev, pch = 16)
  # abline(v = domain[gl_ind], lty = 2)
  # dev.off()

  qsave(bed_sims, file = paste0(train_data_dir, "/bed_sims_", setf, "_", data_date, ".qs"))
  qsave(fric_sims, file = paste0(train_data_dir, "/fric_sims_", setf, "_", data_date, ".qs"))
}

# ## Plot conditional bed simulations
bed_sims_df <- data.frame(
  domain = domain,
  bed1 = bed_sims[, 1],
  bed2 = bed_sims[, 2]
) # ,
# bed_obs = bed_obs_df$bed_elev,
# obs_loc = bed_obs_df$loc) #,

fric_sims_df <- data.frame(
  domain = domain,
  fric1 = fric_sims[, 1],
  fric2 = fric_sims[, 2]
)

bed_sim_plot <- bed_sims_df %>% ggplot() +
  geom_line(aes(x = domain, y = bed1), col = "grey") +
  geom_line(aes(x = domain, y = bed2), col = "grey") +
  geom_point(data = bed_obs_chosen, aes(x = loc, y = bed_elev)) +
  # geom_point(data = bed_obs_val, aes(x = loc, y = bed_elev), col = "red") +
  xlim(c(0, 200e3)) +
  ylim(c(-1500, -500)) +
  theme_bw()

fric_sim_plot <- fric_sims_df %>% ggplot() +
  geom_line(aes(x = domain, y = fric1), col = "grey") +
  geom_line(aes(x = domain, y = fric2), col = "grey") +
  theme_bw()


png(file = paste0("./plots/cnn/input/param_sims_", setf, "_", data_date, ".png"), width = 1000, height = 500)
grid.arrange(bed_sim_plot, fric_sim_plot, nrow = 2)
dev.off()

## Plot parameter simulations
# fric_sims <- qread(file = paste0("./data/training_data/fric_sims_", data_date, ".qs"))
# png(file = paste0("./plots/cnn/input/param_sims_", data_date, ".png"), width = 600, height = 500)
# par(mfrow = c(2,1))
# matplot(bed_sims[, 1:2], type = "l", ylab = "Bedrock Elev", xlab = "Distance along flowline")
# points(bed_obs_chosen$loc, bed_obs_chosen$bed_elev, col = "red", pch = 16)
# matplot(fric_sims[, 1:2], type = "l", ylab = "Friction Coef", xlab = "Distance along flowline")
# dev.off()

# fric_scale <- 1e6 * params$secpera^(1 / params$n)
# param_list <- lapply(1:N, function(r) {
#   list(
#     friction = fric_sims[, r], #* 1e6 * params$secpera^(1 / params$n),
#     bedrock = bed_sims[, r] #+ bed_mean
#   )
# })

## Concatenate the bed simulations and take the mean
# bed_sims_arr <- abind(bed_sim_list, along = 1)
# bed_mean <- colMeans(bed_sims_arr)
