## Generate 1000 simulations 

setwd("/home/babv971/SSA_model/CNN/")

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
regenerate_sims <- F
refit_basis <- F
save_sims <- F
sim_beds <- T

if (sim_beds) {
  train_data_dir <- "./training_data_bed"
} else {
  train_data_dir <- "./training_data"
}

## Presets
data_date <- "20240320" #"20220329" 
N <- 1000 # number of simulations per set
sets <- 1:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 150

# 0. Load ice sheet at steady state
ssa_steady <- readRDS(file = paste0(train_data_dir, "/initial_conds/ssa_steady_20220329.rds", sep = ""))

# 0. Simulate bed observations
set.seed(2024)
ref_bed <- ssa_steady$bedrock
## Randomly select 50 locations between x = 0 and x = 800 km
n_bed_obs <- 50
obs_ind <- sort(sample(length(ref_bed), n_bed_obs))
obs_bed <- ref_bed[obs_ind] + rnorm(n_bed_obs, mean = 0, sd = 20) ## add noise
bed_obs <- list(locations = obs_ind,obs = obs_bed)

years <- 50
domain <- ssa_steady$domain

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1/3)

t1 <- proc.time()
if (regenerate_sims) {

  for (set in sets) {
    cat("Generating set", set, "\n")
    try(generated_data <- run_sims(nsims = N, sim_beds = T, bed_obs = bed_obs)) 

    thickness_velocity_arr_s <- generated_data$thickness_velocity_arr
    friction_arr_s <- generated_data$friction_arr
    gl_arr_s <- generated_data$gl_arr

    if (sim_beds) {
      bed_arr_s <- generated_data$bed_arr
    }

    # ## Basis function representation of the friction coefficients
    # fitted_friction_s <- fit_basis(nbasis = nbasis, domain = domain, friction_arr = friction_arr_s)

    if (save_sims) {
        setf <- formatC(set, width=2, flag="0")
        saveRDS(thickness_velocity_arr_s, file = paste0(train_data_dir, "/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
        saveRDS(friction_arr_s, file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".rds"))
        saveRDS(gl_arr_s, file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".rds"))

        if (sim_beds) {
          saveRDS(bed_arr_s, file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".rds"))

        }
        # saveRDS(fitted_friction_s, file = paste0("./training_data/fitted_friction_", setf, "_", data_date))
    }
    
  }
}


t2 <- proc.time()

## Basis function representation of the friction coefficients
if (refit_basis) {
  for (set in sets) {
    cat("Fitting basis functions for set", set, "\n")
    setf <- formatC(set, width=2, flag="0")
    friction_arr_s <- readRDS(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".rds"))
    bed_arr_s <- readRDS(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".rds"))
    
    ## Should maybe scale the friction values here?
    friction_arr_s <- friction_arr_s/fric_scale
    
    friction_basis <- fit_friction_basis(nbasis = nbasis, domain = domain, sample_arr = friction_arr_s)
    bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, sample_arr = bed_arr_s)

    if (save_sims) {
      saveRDS(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".rds"))
      saveRDS(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".rds"))
      
    }
  }
}

# library(abind)
# f1 <- readRDS(file = paste0("./training_data/friction_basis_", "01", "_", data_date, ".rds"))$fitted_values
# f2 <- readRDS(file = paste0("./training_data/friction_basis_", "02", "_", data_date, ".rds"))$fitted_values
# f12 <- abind(f1, f2, along = 1)
## Now re-scale friction and plot it with the original friction

## Plot some simulations to check 
set <- 1
setf <- formatC(set, width=2, flag="0")
thickness_velocity_arr <- readRDS(file = paste0(train_data_dir, "/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
friction_arr <- readRDS(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".rds"))
bed_arr <- readRDS(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".rds"))
gl_arr <- readRDS(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".rds"))
friction_basis <- readRDS(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".rds"))
bed_basis <- readRDS(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".rds"))

fitted_friction <- friction_basis$fitted_values
fitted_bed <- bed_basis$fitted_values

plots <- list()

nsamples <- 1
sims <- sample(1:N, size = nsamples)

space <- domain/1000
time <- 1:years #1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples*4), nsamples, 4, byrow = T)

gl <- ceiling(gl_arr[1,1]/(domain[length(domain)]/1000)*length(domain))

s <- 1
# for (s in 1:nsamples) {
  sim <- sims[[s]]
  ind <- inds[s, ]
  # print(ind)
  # thickness_velocity_arr <- sim_results[[s]]$thickness_velocity_arr
  thickness <- thickness_velocity_arr[sim,,,1]
  velocity <- thickness_velocity_arr[sim,,,2]
  grid_test$thickness <- as.vector(thickness)
  grid_test$velocity <- as.vector(velocity)

  thickness_plot <- ggplot(grid_test) + 
    geom_tile(aes(space, time, fill = thickness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) + 
    theme_bw() +
    labs(fill="Thickness (m)")
  velocity_plot <- ggplot(grid_test) + 
    geom_tile(aes(space, time, fill = velocity)) +
    theme_bw() +
    scale_fill_distiller(palette = "Reds", direction = 1) + 
    labs(fill=bquote('Velocity (m'~a^-1~')'))
  
  friction_sim <- friction_arr[sim, 1:gl]
  fitted_fric_sim <- fitted_friction[sim, 1:gl]
  df <- data.frame(domain = ssa_steady$domain[1:gl]/1000, friction = friction_sim/fric_scale,
                  fitted_fric = fitted_fric_sim)
  friction_plot <- ggplot(df, aes(x = domain, y = friction)) + geom_line() + 
                    geom_line(aes(x = domain, y = fitted_fric), col = "red") +
                    theme_bw() + 
                    xlab("Domain (km)") + ylab(bquote('Friction (M Pa m'^'-1/3'~'a'^'1/3'~')'))
  
  bed_sim <- bed_arr[sim, ]
  fitted_bed_sim <- fitted_bed[sim, ]
  bed_df <- data.frame(domain = ssa_steady$domain/1000, bed = bed_sim)
  bed_plot <- ggplot(bed_df, aes(x = domain, y = bed)) + geom_line() + 
              geom_line(aes(x = domain, y = fitted_bed_sim), col = "red") +
              theme_bw() + xlab("Domain (km)") + ylab(bquote('Bed (m)'))

  plots[[ind[1]]] <- thickness_plot
  plots[[ind[2]]] <- velocity_plot
  plots[[ind[3]]] <- friction_plot
  plots[[ind[4]]] <- bed_plot
  
# }

png(file = paste0("./plots/sim_beds/simulations_", setf, "_", data_date, ".png"), width = 2000, height = 500)
grid.arrange(grobs = plots, nrow = 1, ncol = 4)
dev.off()

## Plot fitted frictions
nsamples <- 10
par(mfrow = c(nsamples/2, 2))

for (sim in 1:nsamples) {
  plot_domain <- 1:1000
  plot(domain[plot_domain]/1000, friction_arr[sim, plot_domain,,]/fric_scale,
       type = "l", lwd = 1.5, xlab = "Domain (km)", ylab = "Friction (unit)")
  # lines(domain[plot_domain]/1000, lmfit$fitted.values[plot_domain]/fric_scale, col = "seagreen", lwd = 1.5)
  
  lines(domain[plot_domain]/1000, fitted_friction[sim, plot_domain], col = "red", lwd = 1.5)
  # legend("topright", legend = c("global basis", "local basis"), col = c("seagreen", "red"), lty = 1, lwd = 1.5)
  # legend("topright", legend = c("original friction", "local basis rep"), col = c("black", "red"), lty = 1, lwd = 1.5)
  
}

















