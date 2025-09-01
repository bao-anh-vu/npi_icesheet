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

# s <- sample_ind

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

# true_bed <- true_bed[sample_ind, ]
# true_fric <- true_fric[sample_ind, ]
# true_gl <- true_gl[sample_ind, ]

## Choose only the samples that we are interested in
# true_thicknesses_ini <- true_thicknesses[sample_ind, , save_points[1]]
# true_thicknesses_mid <- true_thicknesses[sample_ind, , save_points[2]]
# true_thicknesses_fin <- true_thicknesses[sample_ind, , save_points[3]]
# true_velocities_ini <- true_velocities[sample_ind, , save_points[1]]
# true_velocities_mid <- true_velocities[sample_ind, , save_points[2]]
# true_velocities_fin <- true_velocities[sample_ind, , save_points[3]]

true_thicknesses_ini <- true_thicknesses[, , save_points[1]]
true_thicknesses_mid <- true_thicknesses[, , save_points[2]]
true_thicknesses_fin <- true_thicknesses[, , save_points[3]]
true_velocities_ini <- true_velocities[, , save_points[1]]
true_velocities_mid <- true_velocities[, , save_points[2]]
true_velocities_fin <- true_velocities[, , save_points[3]]


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

## Compute 95% credible intervals for bed and friction
cnn.coverage_bed <- c()
cnn.coverage_fric <- c()
for (r in 1:length(sample_ind)) {
    s <- sample_ind[[r]] # index of sample in the test set (out of 500)
    cnn.bed_q <- apply(bed_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
    cnn.fric_q <- apply(fric_samples_ls[[s]], 1, quantile, probs = c(0.025, 0.975))
    # cnn.bed_lq <- cnn.bed_q[1, ]
    # cnn.bed_uq <- cnn.bed_q[2, ]
    # cnn.fric_lq <- cnn.fric_q[1, ]
    # cnn.fric_uq <- cnn.fric_q[2, ]

    bed_pts_within <- sum(true_bed[s, ] > cnn.bed_q[1, ] & true_bed[s, ] < cnn.bed_q[2, ]) # Check if true bed is within the credible interval
    cnn.coverage_bed[r] <- bed_pts_within / length(true_bed[r, ])

    fric_pts_within <- sum(true_fric[s, ] > cnn.fric_q[1, ] & true_fric[s, ] < cnn.fric_q[2, ]) # Check if true bed is within the credible interval
    cnn.coverage_fric[r] <- fric_pts_within / length(true_fric[r, ])
}

## Plot the sample prediction + credible intervals to check visually
# matplot(t(cnn.bed_q), type = "l", col = "red")
# lines(true_bed[r, ])

## Read CNN-EnKF results
print("Reading posterior samples from CNN-EnKF...")
cnn_enkf.output_dir <- as.list(cnn_enkf.output_dir)
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_thickness_sample", sample_ind[s], "_Ne500_", output_date, ".qs"))
cnn.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
vel_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_velocities_sample", sample_ind[s], "_Ne500_", output_date, ".qs"))
cnn.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
cnn.coverage_thickness <- c()
## Note: cnn.thickness is a list of samples, each sample is a list of ensembles at 3 time points
## So the structure is cnn.thickness[[sample]][[time]]
# cnn.thickness <- lapply(cnn.thickness, function(s) lapply(s, as.matrix))

for (r in 1:length(sample_ind)) {
# for (r in sample_ind) {    
    s <- sample_ind[[r]] # index of sample in the test set (out of 500)
    cat("Processing sample", s, "for N-EnKF... \n")
    cnn.ens_t <- cnn.thickness[[r]][[3]] # Take the last time point
    J <- nrow(cnn.ens_t) # Number of grid points
    cnn.covmat <- 1 / (ncol(cnn.ens_t) - 1) * tcrossprod(cnn.ens_t - rowMeans(cnn.ens_t)) # diag(enkf_covmats[[t]])
    cnn.lower <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(cnn.covmat)[1:J])
    cnn.upper <- rowMeans(cnn.ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(cnn.covmat)[1:J])
    cnn.mean_thickness <- rowMeans(cnn.ens_t)

    cnn.thickness_pts_within <- sum(true_thicknesses_fin[s, ] > cnn.lower & true_thicknesses_fin[s, ] < cnn.upper)
    cnn.coverage_thickness[r] <- cnn.thickness_pts_within / length(true_thicknesses_fin[s, ])

    png(paste0("./plots/coverage/cvg_cnn_thickness_sample", s, "_", output_date, ".png"), width = 800, height = 600)
    plot(cnn.upper, type = "l", col = "red", xlab = "Grid point", ylab = "Thickness (m)")
    lines(cnn.lower, col = "red")
    lines(true_thicknesses_fin[s, ])
    dev.off()
}

cnn.coverage_df <- data.frame(bed = mean(cnn.coverage_bed),
                        fric = mean(cnn.coverage_fric),
                        thickness = mean(cnn.coverage_thickness))
print(cnn.coverage_df)

# cnn.velocity <- lapply(cnn.velocity, function(s) lapply(s, as.matrix))
# cnn.mean_velocity <- lapply(cnn.velocity, function(s) lapply(s, rowMeans))

# cnn.mean_vel_ini <- lapply(cnn.mean_velocity, function(s) s[[1]])
# cnn.mean_vel_mid <- lapply(cnn.mean_velocity, function(s) s[[2]])
# cnn.mean_vel_fin <- lapply(cnn.mean_velocity, function(s) s[[3]])

########################################
##       EnKF-StateAug results        ##
########################################

print("Reading posterior samples from EnKF-StateAug...")
Ne <- 1000 # Ensemble size for Aug-EnKF

aug_enkf.output_dir <- as.list(aug_enkf.output_dir)
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_thickness_sample", sample_ind[s], "_Ne", Ne, "_", output_date, ".qs"))
bed_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_bed_sample", sample_ind[s], "_Ne", Ne, "_", output_date, ".qs"))
fric_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_friction_sample", sample_ind[s], "_Ne", Ne, "_", output_date, ".qs"))
# vel_files <- lapply(1:length(sample_ind), function(s) paste0(aug_enkf.output_dir[[s]], "/enkf_velocities_sample", sample_ind[s], "_Ne", Ne, "_", output_date,  ".qs", sep = ""))

aug_enkf.thickness <- lapply(thickness_files, qread) 
aug_enkf.bed <- lapply(bed_files, qread) 
aug_enkf.fric <- lapply(fric_files, qread) 
# aug_enkf.velocity <- lapply(vel_files, qread) 

# aug_enkf.thickness <- lapply(aug_enkf.thickness, function(s) lapply(s, as.matrix))
# aug_enkf.thickness_fin <- lapply(aug_enkf.thickness, function(s) s[[3]])

aug_enkf.coverage_thickness <- c()
aug_enkf.coverage_bed <- c()
aug_enkf.coverage_fric <- c()
for (r in 1:length(sample_ind)) {
# for (r in sample_ind) {
    s <- sample_ind[[r]] # index of sample in the test set (out of 500)
    cat("Processing sample", s, "for Aug-EnKF... \n")
    aug_enkf.thickness_t <- aug_enkf.thickness[[r]][[3]] # Take the last time point
    aug_enkf.bed_t <- aug_enkf.bed[[r]][[3]]
    aug_enkf.fric_t <- aug_enkf.fric[[r]][[3]]

    aug_enkf.ens_t <- rbind(aug_enkf.thickness_t, aug_enkf.bed_t, aug_enkf.fric_t) # Combine all state variables
    J <- nrow(aug_enkf.thickness_t) # Number of grid points
    aug_enkf.covmat <- 1 / (ncol(aug_enkf.ens_t) - 1) * tcrossprod(aug_enkf.ens_t - rowMeans(aug_enkf.ens_t)) # diag(enkf_covmats[[t]])
    aug_enkf.lower <- rowMeans(aug_enkf.ens_t) + qnorm(0.025) * sqrt(diag(aug_enkf.covmat))
    aug_enkf.upper <- rowMeans(aug_enkf.ens_t) + qnorm(0.975) * sqrt(diag(aug_enkf.covmat))

    aug_enkf.thickness_upper <- aug_enkf.upper[1:J]
    aug_enkf.thickness_lower <- aug_enkf.lower[1:J]
    aug_enkf.bed_upper <- aug_enkf.upper[J + 1:J]
    aug_enkf.bed_lower <- aug_enkf.lower[J + 1:J]
    aug_enkf.fric_upper <- exp(aug_enkf.upper[2*J + 1:J])
    aug_enkf.fric_lower <- exp(aug_enkf.lower[2*J + 1:J])

    aug_enkf.mean_thickness <- rowMeans(aug_enkf.thickness_t)
    aug_enkf.mean_bed <- rowMeans(aug_enkf.bed_t)
    aug_enkf.mean_fric <- rowMeans(exp(aug_enkf.fric_t))

    aug_enkf.thickness_pts_within <- sum(true_thicknesses_fin[s, ] > aug_enkf.thickness_lower & true_thicknesses_fin[s, ] < aug_enkf.thickness_upper)
    aug_enkf.coverage_thickness[r] <- aug_enkf.thickness_pts_within / length(true_thicknesses_fin[s, ])

    aug_enkf.bed_pts_within <- sum(true_bed[s, ] > aug_enkf.bed_lower & true_bed[s, ] < aug_enkf.bed_upper)
    aug_enkf.coverage_bed[r] <- aug_enkf.bed_pts_within / length(true_bed[s, ])
    
    aug_enkf.fric_pts_within <- sum(true_fric[s, ] > aug_enkf.fric_lower & true_fric[s, ] < aug_enkf.fric_upper)
    aug_enkf.coverage_fric[r] <- aug_enkf.fric_pts_within / length(true_fric[s, ])

    ## Plot to check
    png(paste0("./plots/coverage/cvg_aug_enkf_thickness_sample", s, "_Ne", Ne, "_", output_date, ".png"), width = 800, height = 600)
    plot(aug_enkf.thickness_upper, type = "l", col = "red", xlab = "Grid point", ylab = "Thickness (m)")
    lines(aug_enkf.thickness_lower, col = "red")
    lines(true_thicknesses_fin[s, ])
    dev.off()

    png(paste0("./plots/coverage/cvg_aug_enkf_bed_sample", s, "_Ne", Ne, "_", output_date, ".png"), width = 800, height = 600)
    plot(aug_enkf.bed_upper, type = "l", col = "red", xlab = "Grid point", ylab = "Bed (m)")
    lines(aug_enkf.bed_lower, col = "red")
    lines(true_bed[s, ])
    dev.off()

    png(paste0("./plots/coverage/cvg_aug_enkf_fric_sample", s, "_Ne", Ne, "_", output_date, ".png"), width = 800, height = 600)
    plot(aug_enkf.fric_upper, type = "l", col = "red", xlab = "Grid point", ylab = "Friction coeff")
    lines(aug_enkf.fric_lower, col = "red")
    lines(true_fric[s, ])
    dev.off()

}

aug_enkf.coverage_df <- data.frame(bed = mean(aug_enkf.coverage_bed),
                            fric = mean(aug_enkf.coverage_fric),
                            thickness = mean(aug_enkf.coverage_thickness))

## Compare coverage of CNN-EnKF, and Aug-EnKF
coverage_df <- cbind(t(cnn.coverage_df), t(aug_enkf.coverage_df))
colnames(coverage_df) <- c("N-EnKF", "Aug-EnKF")
print(coverage_df)
# aug_enkf.mean_thickness <- lapply(aug_enkf.thickness, function(s) lapply(s, rowMeans))

# aug_enkf.velocity <- lapply(aug_enkf.velocity, function(s) lapply(s, as.matrix))
# aug_enkf.mean_velocity <- lapply(aug_enkf.velocity, function(s) lapply(s, rowMeans))

## Extract mean thicknesses at the beginning, middle, and end of the simulation
# aug_enkf.mean_thickness_ini <- lapply(aug_enkf.mean_thickness, function(s) s[[1]])
# aug_enkf.mean_thickness_mid <- lapply(aug_enkf.mean_thickness, function(s) s[[2]])
# aug_enkf.mean_thickness_fin <- lapply(aug_enkf.mean_thickness, function(s) s[[3]])

# aug_enkf.mean_vel_ini <- lapply(aug_enkf.mean_velocity, function(s) s[[1]])
# aug_enkf.mean_vel_mid <- lapply(aug_enkf.mean_velocity, function(s) s[[2]])
# aug_enkf.mean_vel_fin <- lapply(aug_enkf.mean_velocity, function(s) s[[3]])

## Reorganise these lists so that it has the structure
## cnn.thickness_<time>[grid_pt][sample] 
# aug_enkf.mean_thickness_ini_mat <- t(matrix(unlist(aug_enkf.mean_thickness_ini), nrow = length(cnn.mean_thickness_ini[[1]]), byrow = F))
# aug_enkf.mean_thickness_mid_mat <- t(matrix(unlist(aug_enkf.mean_thickness_mid), nrow = length(cnn.mean_thickness_mid[[1]]), byrow = F))
# aug_enkf.mean_thickness_fin_mat <- t(matrix(unlist(aug_enkf.mean_thickness_fin), nrow = length(cnn.mean_thickness_fin[[1]]), byrow = F))
# aug_enkf.mean_vel_ini_mat <- t(matrix(unlist(aug_enkf.mean_vel_ini), nrow = length(cnn.mean_vel_ini[[1]]), byrow = F))
# aug_enkf.mean_vel_mid_mat <- t(matrix(unlist(aug_enkf.mean_vel_mid), nrow = length(cnn.mean_vel_mid[[1]]), byrow = F))
# aug_enkf.mean_vel_fin_mat <- t(matrix(unlist(aug_enkf.mean_vel_fin), nrow = length(cnn.mean_vel_fin[[1]]), byrow = F))
