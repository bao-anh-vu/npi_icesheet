## Test simulation for a specific bedrock and friction profile

setwd("~/SSA_model/CNN/pilot/simbed")

rm(list = ls())

library("parallel")
library("Matrix")
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library("R.utils")
# library("sp")
library("fields")
# library("tidyr")
library("dplyr")
library("matrixStats") # for the rowMaxs() function
library("mvtnorm")
library("parallel")

library("ggplot2")
# library("plotly")
library("gridExtra")
library(FRK)

# library(gstat)
# library(sp)
# library("mvnfast")
# library("splines")

source("./source/run_sims.R")
source("./source/process_sim_results.R")
source("./source/fit_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
source("./source/create_ref.R")
source("./source/solve_ssa_nl.R")
# source("./source/solve_ssa_nl2.R")

source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/mvnorm_sample.R")
source("./source/get_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")
# source("initialise_ens.R")
# source("ssa_plot_ini_ens.R")
source("./source/propagate.R")
source("./source/obs_operator.R")
# source("run_enkf.R")
# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
# source("initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# ## Seed for generating bed
# ssa_seed <- 123
# set.seed(ssa_seed)

## Some flags
regenerate_sims <- T
refit_basis <- T
save_sims <- T
# sim_beds <- T

# if (sim_beds) {
  train_data_dir <- "~/SSA_model/CNN/pilot/simbed/training_data"
# } else {
#   train_data_dir <- "./training_data"
# }

## Presets
data_date <- "20240320" #"20220329" 
N <- 5000 # number of simulations per set
sets <- 1 #1:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 150

# 0. Load ice sheet at steady state
ssa_steady <- readRDS(file = paste0("~/SSA_model/CNN/training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))

# 0. Simulate bed observations
set.seed(2024)
ref_bed <- ssa_steady$bedrock

## Randomly select 50 locations between x = 0 and x = 800 km
n_bed_obs <- 50
obs_ind <- sort(sample(length(ref_bed), n_bed_obs))
obs_bed <- ref_bed[obs_ind] + rnorm(n_bed_obs, mean = 0, sd = 20) ## add noise
bed_obs <- list(locations = obs_ind,obs = obs_bed)

years <- 20
domain <- ssa_steady$domain

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1/3)

sim_results <- readRDS(file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".rds")) 

## Cores 6 and 2 had errors out of 10 cores, so seems like 20% of the simulations had errors
## Let's see what the errors were
bad_sims <- sim_results$bad_sims
sim_results$errors[[1]]
ind <- bad_sims[1]
bad_bed <- sim_results$params[[ind]]$bedrock
bad_fric <- sim_results$params[[ind]]$friction

good_sims <- sim_results$results[-bad_sims]

## Maybe re-run this simulation manually?
## Process noise parameters
ones <- rep(1, length(domain))
D <- rdist(domain)
l <- 50e3
R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
L <- t(chol(R))
L <- as(L, "dgCMatrix")
process_noise_info <- list(corrmat_chol = L, length_scale = l)

reference <- solve_ssa_nl(
                        domain = ssa_steady$domain,
                        bedrock = bad_bed, 
                        friction_coef = bad_fric,
                        ini_velocity = ssa_steady$current_velocity,
                        ini_thickness = ssa_steady$current_thickness,
                        years = years, steps_per_yr = 52,
                        save_model_output = TRUE,
                        perturb_hardness = TRUE,
                        add_process_noise = T,
                        process_noise_info = process_noise_info
                        )
        
    
            thickness_velocity_obs <- array(
                data = cbind(
                    reference$all_thicknesses[, 2:(years + 1)],
                    reference$all_velocities[, 2:(years + 1)]
                ),
                dim = c(length(domain), years, 2)
            )
            gl <- reference$grounding_line
        # )    

            simulated_data <- list(
                thickness_velocity_arr = thickness_velocity_obs,
                friction_arr = bad_fric,
                bed_arr = bad_bed,
                grounding_line = gl
            )
