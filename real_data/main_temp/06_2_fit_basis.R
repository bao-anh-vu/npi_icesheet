setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library(parallel)
# library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library(R.utils)
# library("sp")
# library(fields)
# library("tidyr")
library(dplyr)
# library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(abind)
library(ggplot2)
library(gridExtra)
library(FRK)
library(qs)

source("./source/fit_basis.R")

## Directory for training data
train_data_dir <- "./data/training_data"
data_date <- "20241111" #"20241103" 

## Flowline data
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Basis function parameters
n_fric_basis <- 120
n_bed_basis <- 150

## Read simulated parameters
sets <- 11#:20

# bed_sim_list <- list()
# fric_sim_list <- list()
# for (i in 1:length(sets)) {
#     setf <- formatC(sets[i], width = 2, flag = "0")
#     bed_sim_list[[i]] <- qread(file = paste0(train_data_dir, "/bed_sims_", setf, "_", data_date, ".qs"))
#     fric_sim_list[[i]] <- qread(file = paste0(train_data_dir, "/fric_sims_", setf, "_", data_date, ".qs"))
# }

# ## Concatenate the bed simulations and take the mean
# bed_sims_arr <- abind(bed_sim_list, along = 1)
# bed_mean <- colMeans(bed_sims_arr)

## Bed prior from GP fit to BedMap data
bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
# L <- t(chol(bed_prior$cov))
bed_mean <- bed_prior$mean
# mean_mat <- matrix(rep(bed_prior$mean, N), nrow = nrow(bed_prior$mean), ncol = N)


## Bedmachine data to compare
bedmachine_data <- qread(paste0("./data/bedmachine/flowline_bedmachine.qs"))
bedmachine <- bedmachine_data$bed_avg

## Then fit basis functions
for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Fitting basis functions for set", set, "\n")
    setf <- formatC(set, width = 2, flag = "0")

    # bed_sims <- bed_sim_list[[i]]
    # fric_sims <- fric_sim_list[[i]]

    # bed_sims <- qread(file = paste0(train_data_dir, "/bed_sims_", setf, "_", data_date, ".qs"))
    # fric_sims <- qread(file = paste0(train_data_dir, "/fric_sims_", setf, "_", data_date, ".qs"))
    bed_sims <- t(qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs")))
    fric_sims <- t(qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs")))


    ## Fit basis to log(friction)
    friction_basis <- fit_friction_basis(
        nbasis = n_fric_basis,
        domain = domain,
        fric_arr = t(fric_sims),
        log_transform = T,
        lengthscale = 4e3
    )

    ## De-mean the bedrock
    ## Replicate the mean bed into a matrix so it's easy to subtract from each bed simulation
    N <- ncol(bed_sims)
    mean_mat <- matrix(rep(bed_mean, N), nrow = J, ncol = N) 

    bed_arr_demean <- bed_sims - mean_mat

    ## Fit basis to de-meaned bedrock
    bed_basis <- fit_bed_basis(
        nbasis = n_bed_basis, domain = domain,
        bed_arr = t(bed_arr_demean),
        lengthscale = 2.5e3
    )
    bed_basis$mean <- bed_mean

    ## Save the fitted basis and the fitted parameters
    qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
    qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
    
}

## Bedmachine data to compare
bedmachine_data <- qread(paste0("./data/bedmachine/flowline_bedmachine.qs"))
bedmachine <- bedmachine_data$bed_avg

## BedMap observations
bed_obs_df <- qread(file = paste0("./data/bedmap/bed_obs_df_all.qs"))

## Plot the fitted friction and bed
set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")
sim <- sample(1:N, 1)
png(file = paste0("./plots/cnn/input/basis_set", setf, "_", data_date, ".png"), width = 800, height = 800, res = 100)
par(mfrow = c(2, 1))
plot_domain <- 1:J # 1000
plot(domain[plot_domain], friction_basis$true_vals[sim, plot_domain], type = "l",
    main = "Simulated friction vs fitted", 
    ylab = "Friction coefficient", 
    xlab = "Distance along flowline (m)")
lines(domain[plot_domain], friction_basis$fitted_values[sim, plot_domain], col = "red")

# png(file = paste0("./plots/cnn/input/bed_basis_", setf, "_", data_date, ".png"), width = 800, height = 400, res = 100)
# plot_domain <- 1:J
plot(domain[plot_domain], bed_basis$true_vals[sim, plot_domain] + bed_mean, type = "l",
    main = "Simulated bed vs BedMachine", ylab = "Elevation (m)", xlab = "Distance along flowline (m)")
lines(domain[plot_domain], bed_basis$fitted_values[sim, plot_domain] + bed_mean, col = "red")
lines(domain[plot_domain], bedmachine[plot_domain], col = "blue", lty = 2)
points(bed_obs_df$loc, bed_obs_df$bed_elev, pch = 16, cex = 0.5, col = "cyan")
legend("topleft", legend = c("Simulated bed", "Fitted bed", "BedMachine"), col = c("black", "red", "blue"), lty = c(1, 1, 2))
dev.off()   

## Plot basis function coefficients
png(file = paste0("./plots/cnn/input/basis_coefs_set", setf, "_", data_date, ".png"), width = 800, height = 800, res = 100)
par(mfrow = c(2, 1))
plot(friction_basis$basis_coefs[sim, ], type = "l", lty = 1,
    main = "Friction basis coefficients",
    ylab = "Coefficient value", xlab = "Simulation index")
plot(bed_basis$basis_coefs[sim, ], type = "l", lty = 1,
    main = "Bed basis coefficients",
    ylab = "Coefficient value", xlab = "Simulation index")
dev.off()