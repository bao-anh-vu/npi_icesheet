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
nbasis <- 150

## Read simulated parameters
sets <- 1:10

bed_sim_list <- list()
fric_sim_list <- list()
for (i in 1:length(sets)) {
    setf <- formatC(sets[i], width = 2, flag = "0")
    bed_sim_list[[i]] <- qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
    fric_sim_list[[i]] <- qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
}

## Concatenate the bed simulations and take the mean
bed_sims_arr <- abind(bed_sim_list, along = 1)
bed_mean <- colMeans(bed_sims_arr)

## Replicate the mean bed into a matrix so it's easy to subtract from each bed simulation
N <- nrow(bed_sim_list[[1]])
mean_mat <- matrix(rep(bed_mean, N), nrow = N, ncol = J, byrow = T) 

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
        log_transform = T,
        lengthscale = 5e3
    )

    ## De-mean the bedrock
    # bed_mean <- colMeans(bed_sims)
    bed_arr_demean <- bed_sims - mean_mat

    ## Fit basis to de-meaned bedrock
    bed_basis <- fit_bed_basis(
        nbasis = nbasis, domain = domain,
        bed_arr = bed_arr_demean,
        lengthscale = 2.5e3
    )
    bed_basis$mean <- bed_mean

    ## Pair the bed-friction observations into a list
    fitted_param_list <- lapply(1:N, function(r) {
        list(
            friction = exp(friction_basis$fitted_values[r, ]),
            bedrock = bed_basis$fitted_values[r, ] + bed_mean
        )
    })

    ## Save the fitted basis and the fitted parameters
    qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
    qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
      
    
}

## Plot the fitted friction and bed
set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")
png(file = paste0("./plots/cnn/input/friction_basis_", setf, "_", data_date, ".png"), width = 800, height = 400)
plot_domain <- 1:J # 1000
plot(domain[plot_domain], friction_basis$true_vals[1, plot_domain], type = "l")
lines(domain[plot_domain], friction_basis$fitted_values[1, plot_domain], col = "red")
dev.off()

png(file = paste0("./plots/cnn/input/bed_basis_", setf, "_", data_date, ".png"))
plot_domain <- 1:J
plot(domain[plot_domain], bed_basis$true_vals[1, plot_domain], type = "l")
lines(domain[plot_domain], bed_basis$fitted_values[1, plot_domain], col = "red")
dev.off()   