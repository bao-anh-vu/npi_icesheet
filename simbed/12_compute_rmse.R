# Calculate RMSE over multiple samples
library(ggplot2)
library(mvtnorm)
library(qs)

setwd("/home/babv971/SSA_model/CNN/simbed/")
rm(list = ls())

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"

use_missing_pattern <- T
use_cov_taper <- F
# use_basis_funs <- T
# plot_ice_thickness <- T
# plot_velocity <- T
# plot_bed <- T
# plot_friction <- T
# plot_gl <- T
save_plots <- T

rmse <- function(estimated, true) {
    stopifnot(length(estimated) == length(true))
    sqrt(mean((estimated - true)^2))
}

## 1. Read samples
sample_ind <- 1:10 #c(1:4, 6:7, 10:15) # test samples to compare
# s <- sample_ind
# set.seed(2024)
# chosen_test_samples <- sample(1:500, 50)
# set.seed(NULL)
# s <- chosen_test_samples[sample_ind] # the actual number of the sample in the test set

years <- 20
save_points <- c(1, floor(years/2) + 1, years+1) #c(1, 11, 21)

## Read test data
sets <- 1:50 # datasets
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

if (use_missing_pattern) {
    data_dir <- paste0("./training_data/", setsf, "/missing")
} else {
    data_dir <- paste0("./training_data/", setsf, "/nonmissing")
}

test_data <- qread(file = paste0(data_dir, "/test_data_", output_date, ".qs"))

true_surface_elevs <- test_data$true_surface_elevs_test
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- exp(test_data$true_fric)
true_gl <- test_data$grounding_line * test_data$sd_gl + test_data$mean_gl

true_bed <- true_bed[sample_ind, ]
true_fric <- true_fric[sample_ind, ]
true_gl <- true_gl[sample_ind, ]

## Choose only the samples that we are interested in
true_thicknesses_ini <- true_thicknesses[sample_ind, , save_points[1]]
true_thicknesses_mid <- true_thicknesses[sample_ind, , save_points[2]]
true_thicknesses_fin <- true_thicknesses[sample_ind, , save_points[3]]
true_velocities_ini <- true_velocities[sample_ind, , save_points[1]]
true_velocities_mid <- true_velocities[sample_ind, , save_points[2]]
true_velocities_fin <- true_velocities[sample_ind, , save_points[3]]

## Read SSA model domain
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
domain <- ssa_steady$domain

## To calculate RMSE for friction, need true gl for each sample
true_gl_pos <- true_gl[, years] * 1000 # Take GL position at last time point for each sample
true_gl_ind <- sapply(1:length(true_gl_pos), function(x) which(domain == true_gl_pos[x]))#floor(true_gl_pos / 800 * 2001) # Convert to grid index

## Read CNN predictions
if (use_missing_pattern) {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/missing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/missing/sample", sample_ind)
    aug_enkf.output_dir <- paste0("./output/stateaug/", setsf, "/missing/sample", sample_ind)    
} else {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing/sample", sample_ind)
    aug_enkf.output_dir <- paste0("./output/stateaug/", setsf, "/nonmissing/sample", sample_ind)    
}

if (use_cov_taper) {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/taper")
    aug_enkf.output_dir <- paste0(aug_enkf.output_dir, "/taper")
} else {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/no_taper")
    aug_enkf.output_dir <- paste0(aug_enkf.output_dir, "/no_taper")
}

pred_fric <- qread(file = paste0(cnn.output_dir, "/pred_fric_", output_date, ".qs"))
pred_bed <- qread(file = paste0(cnn.output_dir, "/pred_bed_", output_date, ".qs"))
pred_gl <- qread(file = paste0(cnn.output_dir, "/pred_gl_", output_date, ".qs"))

print("Reading posterior samples from CNN...")
# fric_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/fric_post_samples_", output_date, ".rds"))
# bed_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/bed_post_samples_", output_date, ".rds"))
# gl_samples_ls <- readRDS(file = paste0(cnn.output_dir, "/gl_post_samples_", output_date, ".rds"))

fric_samples_ls <- qread(file = paste0(cnn.output_dir, "/fric_post_samples_", output_date, ".qs"))
bed_samples_ls <- qread(file = paste0(cnn.output_dir, "/bed_post_samples_", output_date, ".qs"))
gl_samples_ls <- qread(file = paste0(cnn.output_dir, "/gl_post_samples_", output_date, ".qs"))


## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1 / 3)

## Read EnKF-CNN results
print("Reading posterior samples from EnKF-CNN...")
cnn_enkf.output_dir <- as.list(cnn_enkf.output_dir)
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_thickness_sample", s, "_Ne500_", output_date, ".qs"))
cnn.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
vel_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_velocities_sample", s, "_Ne500_", output_date, ".qs"))
cnn.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))

# CNN prediction of bed, friction, GL
cnn.bed <- t(pred_bed[, sample_ind])
cnn.fric <- t(pred_fric[, sample_ind])
cnn.gl <- t(pred_gl[sample_ind, ])

## Note: cnn.thickness is a list of samples, each sample is a list of ensembles at 3 time points
## So the structure is cnn.thickness[[sample]][[time]]
cnn.thickness <- lapply(cnn.thickness, function(s) lapply(s, as.matrix))
cnn.mean_thickness <- lapply(cnn.thickness, function(s) lapply(s, rowMeans))

cnn.velocity <- lapply(cnn.velocity, function(s) lapply(s, as.matrix))
cnn.mean_velocity <- lapply(cnn.velocity, function(s) lapply(s, rowMeans))

## Extract mean thicknesses at the beginning, middle, and end of the simulation
cnn.mean_thickness_ini <- lapply(cnn.mean_thickness, function(s) s[[1]])
cnn.mean_thickness_mid <- lapply(cnn.mean_thickness, function(s) s[[2]])
cnn.mean_thickness_fin <- lapply(cnn.mean_thickness, function(s) s[[3]])

cnn.mean_vel_ini <- lapply(cnn.mean_velocity, function(s) s[[1]])
cnn.mean_vel_mid <- lapply(cnn.mean_velocity, function(s) s[[2]])
cnn.mean_vel_fin <- lapply(cnn.mean_velocity, function(s) s[[3]])

## Reorganise these lists so that it has the structure
## cnn.thickness[[time]][grid_pt][sample] 

cnn.mean_thickness_ini_mat <- t(matrix(unlist(cnn.mean_thickness_ini), nrow = length(cnn.mean_thickness_ini[[1]]), byrow = F))
cnn.mean_thickness_mid_mat <- t(matrix(unlist(cnn.mean_thickness_mid), nrow = length(cnn.mean_thickness_mid[[1]]), byrow = F))
cnn.mean_thickness_fin_mat <- t(matrix(unlist(cnn.mean_thickness_fin), nrow = length(cnn.mean_thickness_fin[[1]]), byrow = F))
cnn.mean_vel_ini_mat <- t(matrix(unlist(cnn.mean_vel_ini), nrow = length(cnn.mean_vel_ini[[1]]), byrow = F))
cnn.mean_vel_mid_mat <- t(matrix(unlist(cnn.mean_vel_mid), nrow = length(cnn.mean_vel_mid[[1]]), byrow = F))
cnn.mean_vel_fin_mat <- t(matrix(unlist(cnn.mean_vel_fin), nrow = length(cnn.mean_vel_fin[[1]]), byrow = F))

## Calculate RMSE (over all samples) for each grid point
cnn.thickness_rmse_ini <- c()
cnn.thickness_rmse_mid <- c()
cnn.thickness_rmse_fin <- c()
cnn.vel_rmse_ini <- c()
cnn.vel_rmse_mid <- c()
cnn.vel_rmse_fin <- c()
for (grid_pt in 1:ncol(cnn.mean_vel_ini_mat)) {
    cnn.thickness_rmse_ini[grid_pt] <- rmse(cnn.mean_thickness_ini_mat[, grid_pt], true_thicknesses_ini[, grid_pt])
    cnn.thickness_rmse_mid[grid_pt] <- rmse(cnn.mean_thickness_mid_mat[, grid_pt], true_thicknesses_mid[, grid_pt])
    cnn.thickness_rmse_fin[grid_pt] <- rmse(cnn.mean_thickness_fin_mat[, grid_pt], true_thicknesses_fin[, grid_pt])
    cnn.vel_rmse_ini[grid_pt] <- rmse(cnn.mean_vel_ini_mat[, grid_pt], true_velocities_ini[, grid_pt])
    cnn.vel_rmse_mid[grid_pt] <- rmse(cnn.mean_vel_mid_mat[, grid_pt], true_velocities_mid[, grid_pt])
    cnn.vel_rmse_fin[grid_pt] <- rmse(cnn.mean_vel_fin_mat[, grid_pt], true_velocities_fin[, grid_pt])
}

## RMSE for bed and friction
cnn.bed_rmse <- c()
cnn.fric_rmse <- c()
for (grid_pt in 1:ncol(cnn.bed)) {
    cnn.bed_rmse[grid_pt] <- rmse(cnn.bed[, grid_pt], true_bed[, grid_pt])
    cnn.fric_rmse[grid_pt] <- rmse(cnn.fric[, grid_pt], true_fric[, grid_pt])
}

## Save plots of RMSE for CNN
if (save_plots) {
    png("./plots/combined/test.png")
    plot(cnn.mean_thickness_fin_mat[1, ], type = "l")
    lines(true_thicknesses[sample_ind[1], , 21], col = "red")
    dev.off()

    png("./plots/combined/rmse_thickness.png", width = 1000, height = 500)
    plot_range <- 1:1000# 2001
    plot(cnn.thickness_rmse_ini[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice thickness")
    lines(cnn.thickness_rmse_mid[plot_range], col = "salmon")
    lines(cnn.thickness_rmse_fin[plot_range], col = "red")
    legend("topright", legend = c("Initial", "Middle", "Final"), col = c("black", "salmon", "red"), lty = 1)
    dev.off()

    png("./plots/combined/rmse_vel.png", width = 1000, height = 500)
    plot_range <- 1:2001
    plot(cnn.vel_rmse_ini[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice velocity")
    lines(cnn.vel_rmse_mid[plot_range], col = "salmon")
    lines(cnn.vel_rmse_fin[plot_range], col = "red")
    legend("topright", legend = c("Initial", "Middle", "Final"), col = c("black", "salmon", "red"), lty = 1)
    dev.off()

    png("./plots/combined/rmse_bed.png", width = 1000, height = 500)
    plot(cnn.bed_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
    dev.off()

    png("./plots/combined/rmse_fric.png", width = 1000, height = 500)
    plot(cnn.fric_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of friction coefficient")
    dev.off()
}

## Read EnKF-StateAug results
print("Reading posterior samples from EnKF-StateAug...")
Ne <- 1000 # Ensemble size for Aug-EnKF

aug_enkf.output_dir <- as.list(aug_enkf.output_dir)
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_thickness_sample", s, "_Ne", Ne, "_", output_date, ".qs"))
aug_enkf.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
vel_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
aug_enkf.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))

aug_enkf.thickness <- lapply(aug_enkf.thickness, function(s) lapply(s, as.matrix))
aug_enkf.mean_thickness <- lapply(aug_enkf.thickness, function(s) lapply(s, rowMeans))

aug_enkf.velocity <- lapply(aug_enkf.velocity, function(s) lapply(s, as.matrix))
aug_enkf.mean_velocity <- lapply(aug_enkf.velocity, function(s) lapply(s, rowMeans))

## Extract mean thicknesses at the beginning, middle, and end of the simulation
aug_enkf.mean_thickness_ini <- lapply(aug_enkf.mean_thickness, function(s) s[[1]])
aug_enkf.mean_thickness_mid <- lapply(aug_enkf.mean_thickness, function(s) s[[2]])
aug_enkf.mean_thickness_fin <- lapply(aug_enkf.mean_thickness, function(s) s[[3]])

aug_enkf.mean_vel_ini <- lapply(aug_enkf.mean_velocity, function(s) s[[1]])
aug_enkf.mean_vel_mid <- lapply(aug_enkf.mean_velocity, function(s) s[[2]])
aug_enkf.mean_vel_fin <- lapply(aug_enkf.mean_velocity, function(s) s[[3]])

## Reorganise these lists so that it has the structure
## cnn.thickness_<time>[grid_pt][sample] 
aug_enkf.mean_thickness_ini_mat <- t(matrix(unlist(aug_enkf.mean_thickness_ini), nrow = length(cnn.mean_thickness_ini[[1]]), byrow = F))
aug_enkf.mean_thickness_mid_mat <- t(matrix(unlist(aug_enkf.mean_thickness_mid), nrow = length(cnn.mean_thickness_mid[[1]]), byrow = F))
aug_enkf.mean_thickness_fin_mat <- t(matrix(unlist(aug_enkf.mean_thickness_fin), nrow = length(cnn.mean_thickness_fin[[1]]), byrow = F))
aug_enkf.mean_vel_ini_mat <- t(matrix(unlist(aug_enkf.mean_vel_ini), nrow = length(cnn.mean_vel_ini[[1]]), byrow = F))
aug_enkf.mean_vel_mid_mat <- t(matrix(unlist(aug_enkf.mean_vel_mid), nrow = length(cnn.mean_vel_mid[[1]]), byrow = F))
aug_enkf.mean_vel_fin_mat <- t(matrix(unlist(aug_enkf.mean_vel_fin), nrow = length(cnn.mean_vel_fin[[1]]), byrow = F))

## Calculate RMSE (over all samples) for each grid point
aug_enkf.thickness_rmse_ini <- c()
aug_enkf.thickness_rmse_mid <- c()
aug_enkf.thickness_rmse_fin <- c()
aug_enkf.vel_rmse_ini <- c()
aug_enkf.vel_rmse_mid <- c()
aug_enkf.vel_rmse_fin <- c()
for (grid_pt in 1:ncol(cnn.mean_vel_ini_mat)) {
    aug_enkf.thickness_rmse_ini[grid_pt] <- rmse(aug_enkf.mean_thickness_ini_mat[, grid_pt], true_thicknesses_ini[, grid_pt])
    aug_enkf.thickness_rmse_mid[grid_pt] <- rmse(aug_enkf.mean_thickness_mid_mat[, grid_pt], true_thicknesses_mid[, grid_pt])
    aug_enkf.thickness_rmse_fin[grid_pt] <- rmse(aug_enkf.mean_thickness_fin_mat[, grid_pt], true_thicknesses_fin[, grid_pt])
    aug_enkf.vel_rmse_ini[grid_pt] <- rmse(aug_enkf.mean_vel_ini_mat[, grid_pt], true_velocities_ini[, grid_pt])
    aug_enkf.vel_rmse_mid[grid_pt] <- rmse(aug_enkf.mean_vel_mid_mat[, grid_pt], true_velocities_mid[, grid_pt])
    aug_enkf.vel_rmse_fin[grid_pt] <- rmse(aug_enkf.mean_vel_fin_mat[, grid_pt], true_velocities_fin[, grid_pt])
}

if (save_plots) {
    png("./plots/combined/compare_rmse_thickness.png", width = 1000, height = 500)
    plot_range <- 1:2001
    plot(aug_enkf.thickness_rmse_fin[plot_range], type = "l", lwd = 2,
            xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice thickness")
    lines(cnn.thickness_rmse_fin[plot_range], col = "salmon", lwd = 2)
    legend("topright", legend = c("Aug-EnKF", "CNN"), col = c("black", "salmon"), lty = 1, lwd = 2)
    dev.off()

    png("./plots/combined/compare_rmse_vel.png", width = 1000, height = 500)
    plot_range <- 1:2001
    plot(aug_enkf.vel_rmse_fin[plot_range], type = "l", lwd = 2,
    xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice velocity")
    lines(cnn.vel_rmse_fin[plot_range], col = "salmon", lwd = 2)
    legend("topright", legend = c("Aug-EnKF", "CNN"), col = c("black", "salmon"), lty = 1, lwd = 2)
    dev.off()

}

## Read bed and friction output from Aug-EnKF
bed_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_bed_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
aug_enkf.bed <- lapply(bed_files, qread) #(file = paste0(aug_enkf.output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
aug_enkf.bed_ini <- lapply(aug_enkf.bed, function(s) s[[1]])
aug_enkf.bed_mid <- lapply(aug_enkf.bed, function(s) s[[2]])
aug_enkf.bed_fin <- lapply(aug_enkf.bed, function(s) s[[3]])
aug_enkf.mean_bed_ini <- lapply(aug_enkf.bed_ini, rowMeans)
aug_enkf.mean_bed_mid <- lapply(aug_enkf.bed_mid, rowMeans)
aug_enkf.mean_bed_fin <- lapply(aug_enkf.bed_fin, rowMeans)

fric_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_friction_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
aug_enkf.fric <- lapply(fric_files, qread) #(file = paste0(aug_enkf.output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
aug_enkf.fric_ini <- lapply(aug_enkf.fric, function(s) s[[1]])
aug_enkf.fric_mid <- lapply(aug_enkf.fric, function(s) s[[2]])
aug_enkf.fric_fin <- lapply(aug_enkf.fric, function(s) s[[3]])
aug_enkf.fric_ini <- lapply(aug_enkf.fric_ini, exp)
aug_enkf.fric_mid <- lapply(aug_enkf.fric_mid, exp)
aug_enkf.fric_fin <- lapply(aug_enkf.fric_fin, exp)
aug_enkf.mean_fric_ini <- lapply(aug_enkf.fric_ini, rowMeans)
aug_enkf.mean_fric_mid <- lapply(aug_enkf.fric_mid, rowMeans)
aug_enkf.mean_fric_fin <- lapply(aug_enkf.fric_fin, rowMeans)

## Matrices containing the mean bed and friction coefficients per sample
## basically each row is one sample
aug_enkf.mean_bed_ini_mat <- t(matrix(unlist(aug_enkf.mean_bed_ini), nrow = length(aug_enkf.mean_bed_ini[[1]]), byrow = F))
aug_enkf.mean_bed_mid_mat <- t(matrix(unlist(aug_enkf.mean_bed_mid), nrow = length(aug_enkf.mean_bed_mid[[1]]), byrow = F))
aug_enkf.mean_bed_fin_mat <- t(matrix(unlist(aug_enkf.mean_bed_fin), nrow = length(aug_enkf.mean_bed_fin[[1]]), byrow = F))

aug_enkf.mean_fric_ini_mat <- t(matrix(unlist(aug_enkf.mean_fric_ini), nrow = length(aug_enkf.mean_fric_ini[[1]]), byrow = F))
aug_enkf.mean_fric_mid_mat <- t(matrix(unlist(aug_enkf.mean_fric_mid), nrow = length(aug_enkf.mean_fric_mid[[1]]), byrow = F))
aug_enkf.mean_fric_fin_mat <- t(matrix(unlist(aug_enkf.mean_fric_fin), nrow = length(aug_enkf.mean_fric_fin[[1]]), byrow = F))

# png("./plots/combined/compare_rmse_bed.png", width = 1000, height = 500)
# plot(aug_enkf.mean_bed_fin_mat[, 1], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
# dev.off()

## RMSE for bed and friction
aug_enkf.bed_rmse <- c()
aug_enkf.fric_rmse <- c()
for (grid_pt in 1:ncol(aug_enkf.mean_bed_fin_mat)) {
    aug_enkf.bed_rmse[grid_pt] <- rmse(aug_enkf.mean_bed_fin_mat[, grid_pt], true_bed[, grid_pt])
    aug_enkf.fric_rmse[grid_pt] <- rmse(aug_enkf.mean_fric_fin_mat[, grid_pt], true_fric[, grid_pt])
}

if (save_plots) {
    png("./plots/combined/compare_rmse_bed.png", width = 1000, height = 500)
    plot(aug_enkf.bed_rmse, type = "l", lwd = 2,
    xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
    lines(cnn.bed_rmse, col = "salmon", lwd = 2)
    legend("topright", legend = c("EnKF", "CNN"), col = c("black", "salmon"), lty = 1, lwd = 2)
    dev.off()

    png("./plots/combined/compare_rmse_fric.png", width = 1000, height = 500)
    plot(aug_enkf.fric_rmse[1:min(true_gl_ind)], type = "l", lwd = 2, 
        xlab = "Grid point", ylab = "RMSE", main = "RMSE of friction coefficient")
    lines(cnn.fric_rmse[1:min(true_gl_ind)], col = "salmon", lwd = 2)
    legend("topright", legend = c("EnKF", "CNN"), col = c("black", "salmon"), lty = 1, lwd = 2)
    dev.off()

    png("./plots/combined/compare_fric.png", width = 1000, height = 500)
    plot(aug_enkf.mean_fric_fin_mat[,1], col = "blue", type = "l")
    lines(cnn.fric[,1], col = "red")
    lines(true_fric[1, ], col = "black")
    dev.off()
}

## Calculate overall RMSE
print("Calculating overall RMSE...")
## Dissect the RMSE in a different way: calculate RMSE per-sample then average over samples
cnn.thickness_rmse_per_sample <- c()
aug_enkf.thickness_rmse_per_sample <- c()
cnn.vel_rmse_per_sample <- c()
aug_enkf.vel_rmse_per_sample <- c()
for (r in 1:length(sample_ind)) {
    cnn.thickness_rmse_per_sample[r] <- rmse(cnn.mean_thickness_fin_mat[r, ], true_thicknesses_fin[r, ])
    aug_enkf.thickness_rmse_per_sample[r] <- rmse(aug_enkf.mean_thickness_ini_mat[r, ], true_thicknesses_ini[r, ])
    cnn.vel_rmse_per_sample[r] <- rmse(cnn.mean_vel_fin_mat[r, ], true_velocities_fin[r, ])
    aug_enkf.vel_rmse_per_sample[r] <- rmse(aug_enkf.mean_vel_ini_mat[r, ], true_velocities_ini[r, ])
}

## Or maybe do RMSE averaged over all samples and all grid points?
cnn.thickness_rmse_ini <- rmse(cnn.mean_thickness_ini_mat, true_thicknesses_ini)
aug_enkf.thickness_rmse_ini <- rmse(aug_enkf.mean_thickness_ini_mat, true_thicknesses_ini)
cnn.vel_rmse_ini <- rmse(cnn.mean_vel_ini_mat, true_velocities_ini)
aug_enkf.vel_rmse_ini <- rmse(aug_enkf.mean_vel_ini_mat, true_velocities_ini)

cnn.thickness_rmse_mid <- rmse(cnn.mean_thickness_mid_mat, true_thicknesses_mid)
aug_enkf.thickness_rmse_mid <- rmse(aug_enkf.mean_thickness_mid_mat, true_thicknesses_mid)
cnn.vel_rmse_mid <- rmse(cnn.mean_vel_mid_mat, true_velocities_mid)
aug_enkf.vel_rmse_mid <- rmse(aug_enkf.mean_vel_mid_mat, true_velocities_mid)

cnn.thickness_rmse_fin <- rmse(cnn.mean_thickness_fin_mat, true_thicknesses_fin)
aug_enkf.thickness_rmse_fin <- rmse(aug_enkf.mean_thickness_fin_mat, true_thicknesses_fin)
cnn.vel_rmse_fin <- rmse(cnn.mean_vel_fin_mat, true_velocities_fin)
aug_enkf.vel_rmse_fin <- rmse(aug_enkf.mean_vel_fin_mat, true_velocities_fin)

cnn.bed_rmse <- rmse(cnn.bed, true_bed)
aug_enkf.bed_rmse_ini <- rmse(aug_enkf.mean_bed_ini_mat, true_bed)
aug_enkf.bed_rmse_mid <- rmse(aug_enkf.mean_bed_mid_mat, true_bed)
aug_enkf.bed_rmse_fin <- rmse(aug_enkf.mean_bed_fin_mat, true_bed)

cnn.fric_rmse_per_sample <- c()
aug_enkf.ini_fric_rmse_per_sample <- c()
aug_enkf.mid_fric_rmse_per_sample <- c()
aug_enkf.fin_fric_rmse_per_sample <- c()
for (r in 1:length(sample_ind)) {
    cnn.fric_rmse_per_sample[r] <- rmse(cnn.fric[r, 1:true_gl_ind[r]], true_fric[r, 1:true_gl_ind[r]])
    aug_enkf.ini_fric_rmse_per_sample[r] <- rmse(aug_enkf.mean_fric_ini_mat[r, 1:true_gl_ind[r]], true_fric[r, 1:true_gl_ind[r]])
    aug_enkf.mid_fric_rmse_per_sample[r] <- rmse(aug_enkf.mean_fric_mid_mat[r, 1:true_gl_ind[r]], true_fric[r, 1:true_gl_ind[r]])
    aug_enkf.fin_fric_rmse_per_sample[r] <- rmse(aug_enkf.mean_fric_fin_mat[r, 1:true_gl_ind[r]], true_fric[r, 1:true_gl_ind[r]])
}

## How to calculate the overall RMSE for friction?
cnn.fric_rmse <- mean(cnn.fric_rmse_per_sample) #rmse(t(cnn.fric), true_fric[s, ])
aug_enkf.fric_rmse_ini <- mean(aug_enkf.ini_fric_rmse_per_sample) #rmse(t(aug_enkf.mean_fric_fin_mat), true_fric[s, ])
aug_enkf.fric_rmse_mid <- mean(aug_enkf.mid_fric_rmse_per_sample) #rmse(t(aug_enkf.mean_fric_fin_mat), true_fric[s, ])
aug_enkf.fric_rmse_fin <- mean(aug_enkf.fin_fric_rmse_per_sample) #rmse(t(aug_enkf.mean_fric_fin_mat), true_fric[s, ])

## RMSE table
rmse_table <- data.frame(
    method = c("CNN", "Aug-EnKF"),
    # thickness_ini = c(cnn.thickness_rmse_ini, aug_enkf.thickness_rmse_ini),
    # thickness_mid = c(cnn.thickness_rmse_mid, aug_enkf.thickness_rmse_mid),
    thickness_fin = c(cnn.thickness_rmse_fin, aug_enkf.thickness_rmse_fin),
    # velocity_ini = c(cnn.vel_rmse_ini, aug_enkf.vel_rmse_ini),
    # velocity_mid = c(cnn.vel_rmse_mid, aug_enkf.vel_rmse_mid),
    # velocity_fin = c(cnn.vel_rmse_fin, aug_enkf.vel_rmse_fin),
    bed = c(cnn.bed_rmse, aug_enkf.bed_rmse_fin),
    friction = c(cnn.fric_rmse, aug_enkf.fric_rmse_fin)
)
print(rmse_table)

## Save RMSE table
write.csv(rmse_table, file = "./output/rmse_table.csv", row.names = F)
