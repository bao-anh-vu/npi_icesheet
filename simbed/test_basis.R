## Test different length scales for the friction basis
setwd("~/SSA_model/CNN/simbed")

rm(list = ls())

library("parallel")
# library("Matrix")
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library("R.utils")
# library("sp")
library("fields")
# library("tidyr")
library("dplyr")
# library("matrixStats") # for the rowMaxs() function
library("mvtnorm")
# library(abind)
library("ggplot2")
# library("plotly")
library("gridExtra")
library("FRK")

source("./source/sim_params.R")
# source("./source/run_sims.R")
source("./source/sim_obs.R")
# source("./source/process_sim_results.R")
source("./source/fit_basis.R")
# source("./source/surface_elev.R")
source("./source/create_params.R")
# source("./source/create_ref.R")
# source("./source/solve_ssa_nl.R")
# source("./source/solve_velocity_azm.R")
# source("./source/solve_thickness.R")
# source("./source/get_surface_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

## Some flagsn,
regenerate_sims <- T
refit_basis <- T
save_sims <- T
log_transform <- T
# sim_beds <- T

train_data_dir <- "./training_data"

## Presets
data_date <- "20240320" # "20220329"
N <- 10 # number of simulations per set
sets <- 51 # 2:5 #1:5h
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 150

# 0. Load ice sheet at steady state
ssa_steady <- readRDS(file = paste0("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))

# 0. Simulate bed observations
set.seed(2024)
ref_bed <- ssa_steady$bedrock
## Randomly select 50 locations between x = 0 and x = 800 km
n_bed_obs <- 50
obs_ind <- sort(sample(length(ref_bed), n_bed_obs))
obs_bed <- ref_bed[obs_ind] + rnorm(n_bed_obs, mean = 0, sd = 20) ## add noise
bed_obs <- list(locations = obs_ind, obs = obs_bed)
rm(.Random.seed, envir = globalenv())

years <- 20
domain <- ssa_steady$domain

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1 / 3)

t1 <- proc.time()

flags <- c()


# for (i in 1:length(sets)) {
#   set <- sets[i]
set <- sets[1]
cat("Generating set", set, "\n")

sim_param_list <- sim_params(
  nsims = N, domain = ssa_steady$domain,
  bed_obs = bed_obs
)

## Then fit basis here
cat("Fitting basis functions for set", set, "\n")
setf <- formatC(set, width = 2, flag = "0")
# friction_arr_s <- readRDS(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".rds"))
# bed_arr_s <- readRDS(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".rds"))

lengthscales <- c(5e3, 10e3, 20e3, 50e3)

friction_ls <- list()
for (i in 1:length(lengthscales)) {
  l <- lengthscales[i]

  ## Fit basis to log(friction)
  friction_ls[[i]] <- fit_friction_basis(
    nbasis = nbasis,
    domain = domain,
    fric_arr = sim_param_list$friction,
    log_transform = log_transform,
    lengthscale = l
  )
}

## Now extract each basis

## Then plot
png(paste0("./plots/basis_fit_log.png"), width = 800, height = 1200)
par(mfrow = c(length(lengthscales), 1))
for (i in 1:length(lengthscales)) {
  friction_basis <- friction_ls[[i]]
  plot_domain <- 1:length(ssa_steady$domain) #1000
  plot(friction_basis$true_vals[1, plot_domain], type = "l",
       xlab = "Grid point", ylab = "Friction coefficient",
       main = paste0("Length scale = ", lengthscales[i] / 1e3, " km"))
  lines(friction_basis$fitted_values[1, plot_domain], col = "red")
}
dev.off()

## Then plot
png(paste0("./plots/basis_fit.png"), width = 800, height = 1200)
par(mfrow = c(length(lengthscales), 1))
for (i in 1:length(lengthscales)) {
  friction_basis <- friction_ls[[i]]
  plot_domain <- 1:length(ssa_steady$domain) #1000
  plot(exp(friction_basis$true_vals[1, plot_domain]), type = "l",
       xlab = "Grid point", ylab = "Friction coefficient",
       main = paste0("Length scale = ", lengthscales[i] / 1e3, " km"))
  lines(exp(friction_basis$fitted_values[1, plot_domain]), col = "red")
}
dev.off()

# Might need different length scales for bedrock

# bed_arr_s <- sim_param_list$bedrock
#     bed_mean <- colMeans(bed_arr_s)
#     mean_mat <- matrix(rep(bed_mean, nrow(bed_arr_s)), nrow = nrow(bed_arr_s), ncol = length(bed_mean), byrow = T)
#     bed_arr_demean <- bed_arr_s - mean_mat

#     ## Fit basis to bedrock
#     bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, bed_arr = bed_arr_demean)
#     # bed_fit <- list(mean = bed_mean, basis = bed_basis)
#     bed_basis$mean <- bed_mean
# }
