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

rmse <- function(estimated, true) {
    stopifnot(length(estimated) == length(true))
    sqrt(mean((estimated - true)^2))
}

## 1. Read samples
sample_ind <- 1:5 # test samples to compare
set.seed(2024)
chosen_test_samples <- sample(1:500, 50)
set.seed(NULL)
s <- chosen_test_samples[sample_ind] # the actual number of the sample in the test set

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

## Read CNN predictions
if (use_missing_pattern) {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/missing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/missing/sample", sample_ind)
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/missing/sample", sample_ind)    
} else {
    cnn.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing")
    cnn_enkf.output_dir <- paste0("./output/posterior/", setsf, "/nonmissing/sample", sample_ind)
    enkfsa.output_dir <- paste0("./output/stateaug/", setsf, "/nonmissing/sample", sample_ind)    
}

if (use_cov_taper) {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/taper")
    enkfsa.output_dir <- paste0(enkfsa.output_dir, "/taper")
} else {
    cnn_enkf.output_dir <- paste0(cnn_enkf.output_dir, "/no_taper")
    enkfsa.output_dir <- paste0(enkfsa.output_dir, "/no_taper")
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
thickness_files <- lapply(sample_ind, function(s) paste0(cnn_enkf.output_dir[s], "/enkf_thickness_sample", s, "_Ne500_", output_date, ".qs"))
cnn.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
vel_files <- lapply(sample_ind, function(s) paste0(cnn_enkf.output_dir[s], "/enkf_velocities_sample", s, "_Ne500_", output_date, ".qs"))
cnn.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))

# CNN prediction of bed, friction, GL
cnn.bed <- pred_bed[, s]
cnn.fric <- pred_fric[, s]
cnn.gl <- pred_gl[s, ]

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

cnn.mean_thickness_ini_mat <- matrix(unlist(cnn.mean_thickness_ini), nrow = length(cnn.mean_thickness_ini[[1]]), byrow = F)
cnn.mean_thickness_mid_mat <- matrix(unlist(cnn.mean_thickness_mid), nrow = length(cnn.mean_thickness_mid[[1]]), byrow = F)
cnn.mean_thickness_fin_mat <- matrix(unlist(cnn.mean_thickness_fin), nrow = length(cnn.mean_thickness_fin[[1]]), byrow = F)
cnn.mean_vel_ini_mat <- matrix(unlist(cnn.mean_vel_ini), nrow = length(cnn.mean_vel_ini[[1]]), byrow = F)
cnn.mean_vel_mid_mat <- matrix(unlist(cnn.mean_vel_mid), nrow = length(cnn.mean_vel_mid[[1]]), byrow = F)
cnn.mean_vel_fin_mat <- matrix(unlist(cnn.mean_vel_fin), nrow = length(cnn.mean_vel_fin[[1]]), byrow = F)

## Calculate RMSE (over all samples) for each grid point
cnn.thickness_rmse_ini <- c()
cnn.thickness_rmse_mid <- c()
cnn.thickness_rmse_fin <- c()
cnn.vel_rmse_ini <- c()
cnn.vel_rmse_mid <- c()
cnn.vel_rmse_fin <- c()
for (grid_pt in 1:nrow(cnn.mean_vel_ini_mat)) {
    cnn.thickness_rmse_ini[grid_pt] <- rmse(cnn.mean_thickness_ini_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[1]])
    cnn.thickness_rmse_mid[grid_pt] <- rmse(cnn.mean_thickness_mid_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[2]])
    cnn.thickness_rmse_fin[grid_pt] <- rmse(cnn.mean_thickness_fin_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[3]])
    cnn.vel_rmse_ini[grid_pt] <- rmse(cnn.mean_vel_ini_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[1]])
    cnn.vel_rmse_mid[grid_pt] <- rmse(cnn.mean_vel_mid_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[2]])
    cnn.vel_rmse_fin[grid_pt] <- rmse(cnn.mean_vel_fin_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[3]])
}

## RMSE for bed and friction
cnn.bed_rmse <- c()
cnn.fric_rmse <- c()
for (grid_pt in 1:nrow(cnn.bed)) {
    cnn.bed_rmse[grid_pt] <- rmse(cnn.bed[grid_pt, ], true_bed[s, grid_pt])
    cnn.fric_rmse[grid_pt] <- rmse(cnn.fric[grid_pt, ], true_fric[s, grid_pt])
}

## Then calculate RMSE over all samples for each time point and each grid point

png("./plots/temp/test.png")
plot(cnn.mean_thickness_fin_mat[, 1], type = "l")
lines(true_thicknesses[s[1], , 21], col = "red")
dev.off()

png("./plots/temp/rmse_thickness.png", width = 1000, height = 500)
plot_range <- 1:1000# 2001
plot(cnn.thickness_rmse_ini[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice thickness")
lines(cnn.thickness_rmse_mid[plot_range], col = "salmon")
lines(cnn.thickness_rmse_fin[plot_range], col = "red")
legend("topright", legend = c("Initial", "Middle", "Final"), col = c("black", "salmon", "red"), lty = 1)
dev.off()

png("./plots/temp/rmse_vel.png", width = 1000, height = 500)
plot_range <- 1:2001
plot(cnn.vel_rmse_ini[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice velocity")
lines(cnn.vel_rmse_mid[plot_range], col = "salmon")
lines(cnn.vel_rmse_fin[plot_range], col = "red")
legend("topright", legend = c("Initial", "Middle", "Final"), col = c("black", "salmon", "red"), lty = 1)
dev.off()

png("./plots/temp/rmse_bed.png", width = 1000, height = 500)
plot(cnn.bed_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
dev.off()

png("./plots/temp/rmse_fric.png", width = 1000, height = 500)
plot(cnn.fric_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of friction coefficient")
dev.off()
# cnn.bed_q <- apply(bed_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
# cnn.fric_q <- apply(fric_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
# cnn.bed_lq <- cnn.bed_q[1, ]
# cnn.bed_uq <- cnn.bed_q[2, ]
# cnn.fric_lq <- cnn.fric_q[1, ]
# cnn.fric_uq <- cnn.fric_q[2, ]

## Read EnKF-StateAug results
print("Reading posterior samples from EnKF-StateAug...")
Ne <- 1000 # Ensemble size for EnKFSA

enkfsa.output_dir <- as.list(enkfsa.output_dir)
thickness_files <- lapply(sample_ind, function(s) paste0(enkfsa.output_dir[s], "/enkf_thickness_sample", s, "_Ne", Ne, "_", output_date, ".qs"))
enkfsa.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
vel_files <- lapply(sample_ind, function(s) paste0(enkfsa.output_dir[s], "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))



enkfsa.thickness <- lapply(enkfsa.thickness, function(s) lapply(s, as.matrix))
enkfsa.mean_thickness <- lapply(enkfsa.thickness, function(s) lapply(s, rowMeans))

enkfsa.velocity <- lapply(enkfsa.velocity, function(s) lapply(s, as.matrix))
enkfsa.mean_velocity <- lapply(enkfsa.velocity, function(s) lapply(s, rowMeans))

## Extract mean thicknesses at the beginning, middle, and end of the simulation
enkfsa.mean_thickness_ini <- lapply(enkfsa.mean_thickness, function(s) s[[1]])
enkfsa.mean_thickness_mid <- lapply(enkfsa.mean_thickness, function(s) s[[2]])
enkfsa.mean_thickness_fin <- lapply(enkfsa.mean_thickness, function(s) s[[3]])

enkfsa.mean_vel_ini <- lapply(enkfsa.mean_velocity, function(s) s[[1]])
enkfsa.mean_vel_mid <- lapply(enkfsa.mean_velocity, function(s) s[[2]])
enkfsa.mean_vel_fin <- lapply(enkfsa.mean_velocity, function(s) s[[3]])

## Reorganise these lists so that it has the structure
## cnn.thickness[[time]][grid_pt][sample] 

enkfsa.mean_thickness_ini_mat <- matrix(unlist(enkfsa.mean_thickness_ini), nrow = length(cnn.mean_thickness_ini[[1]]), byrow = F)
enkfsa.mean_thickness_mid_mat <- matrix(unlist(enkfsa.mean_thickness_mid), nrow = length(cnn.mean_thickness_mid[[1]]), byrow = F)
enkfsa.mean_thickness_fin_mat <- matrix(unlist(enkfsa.mean_thickness_fin), nrow = length(cnn.mean_thickness_fin[[1]]), byrow = F)
enkfsa.mean_vel_ini_mat <- matrix(unlist(enkfsa.mean_vel_ini), nrow = length(cnn.mean_vel_ini[[1]]), byrow = F)
enkfsa.mean_vel_mid_mat <- matrix(unlist(enkfsa.mean_vel_mid), nrow = length(cnn.mean_vel_mid[[1]]), byrow = F)
enkfsa.mean_vel_fin_mat <- matrix(unlist(enkfsa.mean_vel_fin), nrow = length(cnn.mean_vel_fin[[1]]), byrow = F)

## Calculate RMSE (over all samples) for each grid point
enkfsa.thickness_rmse_ini <- c()
enkfsa.thickness_rmse_mid <- c()
enkfsa.thickness_rmse_fin <- c()
enkfsa.vel_rmse_ini <- c()
enkfsa.vel_rmse_mid <- c()
enkfsa.vel_rmse_fin <- c()
for (grid_pt in 1:nrow(cnn.mean_vel_ini_mat)) {
    enkfsa.thickness_rmse_ini[grid_pt] <- rmse(enkfsa.mean_thickness_ini_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[1]])
    enkfsa.thickness_rmse_mid[grid_pt] <- rmse(enkfsa.mean_thickness_mid_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[2]])
    enkfsa.thickness_rmse_fin[grid_pt] <- rmse(enkfsa.mean_thickness_fin_mat[grid_pt, ], true_thicknesses[s, grid_pt, save_points[3]])
    enkfsa.vel_rmse_ini[grid_pt] <- rmse(enkfsa.mean_vel_ini_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[1]])
    enkfsa.vel_rmse_mid[grid_pt] <- rmse(enkfsa.mean_vel_mid_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[2]])
    enkfsa.vel_rmse_fin[grid_pt] <- rmse(enkfsa.mean_vel_fin_mat[grid_pt, ], true_velocities[s, grid_pt, save_points[3]])
}

png("./plots/temp/compare_rmse_thickness.png", width = 1000, height = 500)
plot_range <- 1:1000# 2001
plot(enkfsa.thickness_rmse_fin[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice thickness")
lines(cnn.thickness_rmse_fin[plot_range], col = "salmon")
legend("topright", legend = c("EnKFSA", "CNN"), col = c("black", "salmon"), lty = 1)
dev.off()

png("./plots/temp/compare_rmse_vel.png", width = 1000, height = 500)
plot_range <- 1:2001
plot(enkfsa.vel_rmse_fin[plot_range], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of ice velocity")
lines(cnn.vel_rmse_fin[plot_range], col = "salmon")
legend("topright", legend = c("EnKFSA", "CNN"), col = c("black", "salmon"), lty = 1)
dev.off()


bed_files <- lapply(sample_ind, function(s) paste0(enkfsa.output_dir[s], "/enkf_bed_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.bed <- lapply(bed_files, qread) #(file = paste0(enkfsa.output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.bed_fin <- lapply(enkfsa.bed, function(s) s[[3]])
enkfsa.mean_bed_fin <- lapply(enkfsa.bed_fin, rowMeans)

fric_files <- lapply(sample_ind, function(s) paste0(enkfsa.output_dir[s], "/enkf_friction_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.fric <- lapply(fric_files, qread) #(file = paste0(enkfsa.output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.fric_fin <- lapply(enkfsa.fric, function(s) s[[3]])
enkfsa.fric_fin <- lapply(enkfsa.fric_fin, exp)
enkfsa.mean_fric_fin <- lapply(enkfsa.fric_fin, rowMeans)

enkfsa.mean_bed_fin_mat <- matrix(unlist(enkfsa.mean_bed_fin), nrow = length(enkfsa.mean_bed_fin[[1]]), byrow = F)
enkfsa.mean_fric_fin_mat <- matrix(unlist(enkfsa.mean_fric_fin), nrow = length(enkfsa.mean_fric_fin[[1]]), byrow = F)

# png("./plots/temp/compare_rmse_bed.png", width = 1000, height = 500)
# plot(enkfsa.mean_bed_fin_mat[, 1], type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
# dev.off()

## RMSE for bed and friction
enkfsa.bed_rmse <- c()
enkfsa.fric_rmse <- c()
for (grid_pt in 1:nrow(enkfsa.mean_bed_fin_mat)) {
    enkfsa.bed_rmse[grid_pt] <- rmse(enkfsa.mean_bed_fin_mat[grid_pt, ], true_bed[s, grid_pt])
    enkfsa.fric_rmse[grid_pt] <- rmse(enkfsa.mean_fric_fin_mat[grid_pt, ], true_fric[s, grid_pt])
}

png("./plots/temp/compare_rmse_bed.png", width = 1000, height = 500)
plot(enkfsa.bed_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of bed elevation")
lines(cnn.bed_rmse, col = "salmon")
legend("topright", legend = c("EnKF", "CNN"), col = c("black", "salmon"), lty = 1)
dev.off()

png("./plots/temp/compare_rmse_fric.png", width = 1000, height = 500)
plot(enkfsa.fric_rmse, type = "l", xlab = "Grid point", ylab = "RMSE", main = "RMSE of friction coefficient")
lines(cnn.fric_rmse, col = "salmon")
legend("topright", legend = c("EnKF", "CNN"), col = c("black", "salmon"), lty = 1)
dev.off()

png("./plots/temp/compare_fric.png", width = 1000, height = 500)
plot(enkfsa.mean_fric_fin_mat[,1], col = "blue", type = "l")
lines(cnn.fric[,1], col = "red")
lines(true_fric[s[1], ], col = "black")
dev.off()

# Note: modify plot for friction so that RMSE is calculated only up to the grounding line