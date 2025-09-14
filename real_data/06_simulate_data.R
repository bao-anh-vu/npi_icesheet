# NOTE: need to either smooth out the bed simulations or 
# find a way to modify the parameters in the conditional sim so that the bed isn't so noisy

## Generate 1000 simulations

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library(parallel)
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
library(abind)
library(ggplot2)
library(gridExtra)
library(FRK)
library(qs)

source("./source/sim_params.R")
# source("./source/run_sims.R")
source("./source/sim_obs.R")
source("./source/cond_sim_gp.R")
source("./source/get_ini_thickness.R")
source("./source/process_sim_results.R")
source("./source/fit_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
# source("./source/create_ref.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/get_surface_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

## Some flags
regenerate_sims <- T
refit_basis <- T
save_sims <- T
log_transform <- T
use_basal_melt_data <- F

## Directory for training data
train_data_dir <- "./data/training_data"

## Presets
data_date <- "20241111" #"20241103" 
N <- 1000 # number of simulations per set
# set <- 1 #commandArgs(trailingOnly = TRUE)
sets <- 1:10 #50 #:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])
warmup <- 0
years <- 10 + warmup
nbasis <- 100 

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

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# 0. Load surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# 0. Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))

## Scaling units for friction coefficients
fric_scale <- 1e6 * params$secpera^(1 / params$n)

## SMB data
smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)
params$as <- smb_avg # surface accumulation rate (m/s)

# ## Basal melt data (only available on shelves???)
# if (use_basal_melt_data) {
#   melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
#   # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
#   avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
#   melt_nonmissing <- which(!is.na(avg_melt_rate))
#   avg_melt_rate[1:(melt_nonmissing[1]-1)] <- -2 #seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
#   avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) #mean(avg_melt_rate[melt_nonmissing])
#   avg_melt_rate <- - avg_melt_rate # inverting this as eventually smb is calculated as smb - melt

# } else {
#   avg_melt_rate <- rep(1, length(domain))
# }
avg_melt_rate <- rep(0, length(domain)) # assume no melt for now
params$ab <- avg_melt_rate # basal melt rate (m/s)
# plot(avg_melt_rate, type = "l")

t1 <- proc.time()
if (regenerate_sims) {

  flags <- c()

  bed_sim_list <- list()
  fric_sim_list <- list()

  for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Generating set", set, "\n")
    
    ## Simulate bed and friction
    sim_param_output <- sim_params(
      nsims = N, domain = ssa_steady$domain,
      bed_obs = bed_obs_df[bed_obs_df$chosen == 1, ]
    )

    sim_param_list <- sim_param_output$sim_param_list
    bed_sims <- sim_param_list$bedrock
    fric_sims <- sim_param_list$friction

    bed_sim_list[[i]] <- bed_sims
    fric_sim_list[[i]] <- fric_sims

    ## Plot conditional bed simulations
    bed_sims_df <- data.frame(domain = domain, 
                              bed1 = bed_sims[1, ], 
                              bed2 = bed_sims[2, ],
                              bed_obs = bed_obs_df$bed_elev, 
                              obs_loc = bed_obs_df$loc) #,

    png(file = paste0("./plots/bed/bed_sims_", setf, "_", data_date, ".png"), width = 1000, height = 500)
    bed_sim_plot <- bed_sims_df %>% ggplot() + 
      geom_line(aes(x = domain, y = bed1), col = "grey") +
      geom_line(aes(x = domain, y = bed2), col = "grey") +
      geom_point(data = bed_obs_df, aes(x = loc, y = bed_elev)) +
      theme_bw()
    print(bed_sim_plot)
    dev.off()

    if (save_sims) {
      qsave(bed_sims, file = paste0(train_data_dir, "/bed_sims_", data_date, ".qs"))
      qsave(fric_sims, file = paste0(train_data_dir, "/fric_sims_", data_date, ".qs"))
    }
  }

## Plot parameter simulations
# fric_sims <- qread(file = paste0("./data/training_data/fric_sims_", data_date, ".qs"))
png(file = paste0("./plots/temp/param_sims_", data_date, ".png"), width = 600, height = 500)
par(mfrow = c(2,1))
matplot(t(bed_sims[1:3, ]), type = "l", ylab = "Bedrock Elev", xlab = "Distance along flowline")
matplot(t(fric_sims[1:3, ]), type = "l", ylab = "Friction Coef", xlab = "Distance along flowline")
dev.off()

  ## Concatenate the bed simulations and take the mean
  bed_sims_arr <- abind(bed_sim_list, along = 1)
  bed_mean <- colMeans(bed_sims_arr)
    
  ## Then fit basis functions
  for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Fitting basis functions for set", set, "\n")
    setf <- formatC(set, width = 2, flag = "0")

    bed_sims <- bed_sim_list[[i]] 
    fric_sims <- fric_sim_list[[i]] 

    ## Fit basis to log(friction)
    friction_basis <- fit_friction_basis(
      nbasis = nbasis,
      domain = domain,
      fric_arr = fric_sims,
      log_transform = log_transform,
      lengthscale = 10e3
    )

    ## De-mean the bedrock
    # bed_mean <- colMeans(bed_sims)
    mean_mat <- matrix(rep(bed_mean, nrow(bed_sims)), nrow = nrow(bed_sims), ncol = length(bed_mean), byrow = T)
    bed_arr_demean <- bed_sims - mean_mat

    ## Fit basis to de-meaned bedrock
    bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, 
                              bed_arr = bed_arr_demean,
                              lengthscale = 2.5e3) 
    bed_basis$mean <- bed_mean

    png(file = paste0("./plots/friction/friction_basis_", setf, "_", data_date, ".png"))
    plot_domain <- 1:J #1000
    plot(domain[plot_domain], friction_basis$true_vals[2, plot_domain], type = "l")
    lines(domain[plot_domain], friction_basis$fitted_values[2, plot_domain], col  = "red")
    dev.off()

    png(file = paste0("./plots/bed/bed_basis_", setf, "_", data_date, ".png"))
    plot_domain <- 1:J
    plot(domain[plot_domain], bed_basis$true_vals[1, plot_domain], type = "l")
    lines(domain[plot_domain], bed_basis$fitted_values[1, plot_domain], col  = "red")
    dev.off()

    ## Pair the bed-friction observations into a list
    fitted_param_list <- lapply(1:N, function(r) {
      list(
        friction = friction_basis$fitted_values[r, ],
        bedrock = bed_basis$fitted_values[r, ] + bed_mean
      )
    })

    ## Plot the fitted friction and bed


    ## Generate observations based on the simulated bed and friction
    test <- try(
      sim_results <- sim_obs(
        param_list = fitted_param_list,
        domain = domain,
        phys_params = params,
        years = years, # sim_beds = T,
        warmup = warmup,
        ini_thickness = ssa_steady$current_thickness,
        ini_velocity = ssa_steady$current_velocity,
        smb = smb_avg,
        basal_melt = avg_melt_rate,
        log_transform = log_transform
      )
    )

    # calc <- if(length(sim_results$errors) > 0) {
    #   flags[i] <- 1
    #   next
    # } else {
    #   flags[i] <- 0
    # }

    if (save_sims) {
      qsave(sim_results, file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))
    }

    ## Need to get rid of the simulations that failed here
    bad_sims <- sim_results$bad_sims
    if (length(bad_sims) > 0) {
      cat("Some simulations failed in set", set, "\n")
      good_sims <- sim_results$results[-bad_sims]

      ## Delete corresponding bad simulations from the fitted basis
      friction_basis$basis_coefs <- friction_basis$basis_coefs[-bad_sims, ]
      friction_basis$true_vals <- friction_basis$true_vals[-bad_sims, ]
      friction_basis$fitted_values <- friction_basis$fitted_values[-bad_sims, ]

      bed_basis$basis_coefs <- bed_basis$basis_coefs[-bad_sims, ]
      bed_basis$fitted_values <- bed_basis$fitted_values[-bad_sims, ]
      bed_basis$true_vals <- bed_basis$true_vals[-bad_sims, ]

    } else {
      good_sims <- sim_results$results
    }

    generated_data <- process_sim_results(sims = good_sims)

    # thickness_velocity_arr_s <- generated_data$thickness_velocity_arr
    surface_obs_arr_s <- generated_data$surface_obs_arr
    friction_arr_s <- generated_data$friction_arr
    gl_arr_s <- generated_data$gl_arr
    bed_arr_s <- generated_data$bed_arr

    ## Should scale the friction values here
    # friction_arr_s <- friction_arr_s/fric_scale

    true_surface_elevs <- generated_data$true_surface_elevs
    true_thicknesses <- generated_data$true_thicknesses
    true_velocities <- generated_data$true_velocities

    if (save_sims) {
      setf <- formatC(set, width = 2, flag = "0")
      qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
      qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
      qsave(surface_obs_arr_s, file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
      qsave(friction_arr_s, file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
      qsave(gl_arr_s, file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))
      qsave(bed_arr_s, file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
      qsave(true_surface_elevs, file = paste0(train_data_dir, "/true_surface_elevs_", setf, "_", data_date, ".qs"))
      qsave(true_thicknesses, file = paste0(train_data_dir, "/true_thicknesses_", setf, "_", data_date, ".qs"))
      qsave(true_velocities, file = paste0(train_data_dir, "/true_velocities_", setf, "_", data_date, ".qs"))
    }
  }
} else {

  # flags <- c()
  # for (s in 1:length(sets)) {
  #   set <- sets[s]
  #   cat("Checking set", set, "\n")
  #   setf <- formatC(set, width = 2, flag = "0")
  #   sim_results <- qread(file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))

  #   if (length(sim_results$errors) > 0) {
  #     flags[s] <- 1
  #   } else {
  #     flags[s] <- 0
  #   }
  # }

# bad_sets <- sets[which(flags == 1)]
# sink("bad_sets.txt")
# print(bad_sets)
# sink()
}

t2 <- proc.time()
print(t2 - t1)
# true_vels <- qread(file = paste0("data/training_data/true_velocities_02_", data_date, ".qs"))

########################################
##          Plot simulations          ##
########################################

set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")

### True thickness, friction, bed, grounding line
surface_obs_arr <- qread(file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))

## Fitted friction and bed
friction_basis <- qread(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
friction_arr <- friction_basis$true_vals
fitted_friction <- friction_basis$fitted_values

bed_basis <- qread(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
bed_arr <- bed_basis$true_vals #+ bed_mean_mat
bed_mean <- bed_basis$mean
fitted_bed <- bed_basis$fitted_values

# fitted_friction <- fitted_friction[-bad_sims, ]
# fitted_bed <- fitted_bed[-bad_sims, ]

## Check if the generated data is close to real data here
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs")
vel_mat2 <- qread(file = "./data/velocity/all_velocity_arr.qs")
vel_mat <- qread(file = "./data/velocity/vel_smoothed.qs")

sim <- 2
png("./plots/temp/sim_vs_real.png", width = 500, height = 600)

par(mfrow = c(2, 1))

matplot(surface_obs_arr[sim,,,1], type = "l", col = "grey", 
        xlab = "Grid point", ylab = "Surface elevation (m)")
matlines(surf_elev_mat, col = "salmon", lty = 2)
lines(ssa_steady$current_top_surface, col = "blue")

matplot(vel_mat, ylim = c(0, 6000), type = "l", col = "salmon",
        xlab = "Grid point", ylab = "Velocity (m/a)")
# matlines(vel_mat2, col = "black", lty = 2)
matlines(surface_obs_arr[sim,,,2], col = "grey", lty = 2)
lines(ssa_steady$current_velocity, col = "blue")

dev.off()

## Hovmoller plots
plots <- list()

nsamples <- 4
sims <- sample(1:dim(surface_obs_arr)[1], size = nsamples)

space <- domain / 1000
time <- 1:(years + 1) # 1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples * 4), nsamples, 4, byrow = T)

gl <- ceiling(gl_arr[1, 1] / (domain[length(domain)] / 1000) * length(domain))

# s <- 1
for (s in 1:nsamples) {
  sim <- sims[[s]]
  ind <- inds[s, ]
  surface_elev <- surface_obs_arr[sim, , , 1]
  velocity <- surface_obs_arr[sim, , , 2]
  grid_test$surface_elev <- as.vector(surface_elev)
  grid_test$velocity <- as.vector(velocity)

  surface_elev_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = surface_elev)) +
    scale_y_reverse() +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    # labs(fill="Thickness (m)")
    labs(fill = "Surface elev. (m)")

  velocity_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = velocity)) +
    scale_y_reverse() +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")"))


  if (log_transform) {
    fitted_fric_sim <- exp(fitted_friction[sim, 1:gl])
    friction_sim <- exp(friction_arr[sim, 1:gl])
  } else {
    fitted_fric_sim <- fitted_friction[sim, 1:gl]
    friction_sim <- friction_arr[sim, 1:gl]
  }

  df <- data.frame(
    domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim,
    fitted_fric = fitted_fric_sim
  )
  friction_plot <- ggplot(df, aes(x = domain, y = friction)) +
    geom_line() +
    geom_line(aes(x = domain, y = fitted_fric), col = "red") +
    theme_bw() +
    xlab("Domain (km)") +
    ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")"))

  bed_sim <- bed_arr[sim, ] + bed_mean
  fitted_bed_sim <- fitted_bed[sim, ] + bed_mean
  bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim, fitted_bed = fitted_bed_sim)
  
  bed_plot <- ggplot(bed_df, aes(x = domain, y = bed)) +
    geom_line() +
    geom_line(aes(x = domain, y = fitted_bed), col = "red") +
    theme_bw() +
    xlab("Domain (km)") +
    ylab("Bed (m)")

  # plots[[ind[1]]] <- thickness_plot
  plots[[ind[1]]] <- surface_elev_plot
  plots[[ind[2]]] <- velocity_plot
  plots[[ind[3]]] <- friction_plot
  plots[[ind[4]]] <- bed_plot

}

png(file = paste0("./plots/simulations_", setf, "_", data_date, ".png"), 
          width = 2400, height = 400 * nsamples)
grid.arrange(grobs = plots, nrow = nsamples, ncol = 4)
dev.off()

# png(file = paste0("./plots/cnn/test.png"))
# matplot(surface_obs_arr[2,,,1], type = "l", col = "grey")
# matlines(surf_elev_mat, col = "red")
# dev.off()


