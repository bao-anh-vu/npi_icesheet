# Calculate CRPS over multiple samples
library(ggplot2)
library(mvtnorm)
library(qs)

setwd("/home/babv971/SSA_model/CNN/simbed/")
rm(list = ls())

recalculate_crps <- T

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"

use_missing_pattern <- T
use_cov_taper <- F
save_plots <- T

rmse <- function(estimated, true) {
    stopifnot(length(estimated) == length(true))
    sqrt(mean((estimated - true)^2))
}

## 1. Read samples
sample_ind <- 1:10 #c(1:2, 6:7, 10:15) # test samples to compare
s <- sample_ind
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

true_bed <- true_bed[s, ]
true_fric <- true_fric[s, ]
true_gl <- true_gl[s, ]

## Choose only the samples that we are interested in
true_thicknesses_ini <- true_thicknesses[s, , save_points[1]]
true_thicknesses_mid <- true_thicknesses[s, , save_points[2]]
true_thicknesses_fin <- true_thicknesses[s, , save_points[3]]
true_velocities_ini <- true_velocities[s, , save_points[1]]
true_velocities_mid <- true_velocities[s, , save_points[2]]
true_velocities_fin <- true_velocities[s, , save_points[3]]

## Read SSA model domain
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
domain <- ssa_steady$domain

## To calculate RMSE for friction, need true gl for each sample
# true_gl_pos <- true_gl[years] * 1000 # Take GL position at last time point for each sample
# true_gl_ind <- sapply(1:length(true_gl_pos), function(x) which(domain == true_gl_pos[x]))#floor(true_gl_pos / 800 * 2001) # Convert to grid index

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

## True state
true_thicknesses_fin <- true_thicknesses[s, , save_points[3]]

## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1 / 3)

## Read EnKF-CNN results (for state inference)
print("Reading posterior samples from EnKF-CNN...")
cnn_enkf.output_dir <- as.list(cnn_enkf.output_dir)

## Note: cnn.thickness is a list of samples, each sample is a list of ensembles at 3 time points
## So the structure is cnn.thickness[[sample]][[time]]
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_thickness_sample", sample_ind[s], "_Ne500_", output_date, ".qs"))
cnn.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
cnn.thickness <- lapply(cnn.thickness, function(s) lapply(s, as.matrix))
cnn.thickness_fin <- lapply(cnn.thickness, function(s) s[[3]]) ## extract the final time point

vel_files <- lapply(1:length(sample_ind), function(s) paste0(cnn_enkf.output_dir[[s]], "/enkf_velocities_sample", sample_ind[s], "_Ne500_", output_date, ".qs"))
cnn.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
cnn.velocity <- lapply(cnn.velocity, function(s) lapply(s, as.matrix))

# CNN prediction of bed, friction, GL
cnn.bed_samples <- lapply(s, function(ind) bed_samples_ls[[ind]])
cnn.fric_samples <- lapply(s, function(ind) fric_samples_ls[[ind]]) 
cnn.gl_samples <- lapply(s, function(ind) gl_samples_ls[[ind]])

## Read EnKF-StateAug results
print("Reading posterior samples from EnKF-StateAug...")
Ne <- 1000 # Ensemble size for EnKFSA

enkfsa.output_dir <- as.list(enkfsa.output_dir)
thickness_files <- lapply(1:length(sample_ind), function(s) paste0(enkfsa.output_dir[[s]], "/enkf_thickness_sample", sample_ind[s], "_Ne", Ne, "_", output_date, ".qs"))
enkfsa.thickness <- lapply(thickness_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_thickness_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
enkfsa.thickness <- lapply(enkfsa.thickness, function(s) lapply(s, as.matrix))
enkfsa.thickness_fin <- lapply(enkfsa.thickness, function(s) s[[3]])

# vel_files <- lapply(sample_ind, function(s) paste0(enkfsa.output_dir[s], "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
# enkfsa.velocity <- lapply(vel_files, qread) #(file = paste0(cnn_enkf.output_dir, "/enkf_velocities_sample", sample_ind, "_Ne500_", output_date, ".qs", sep = ""))
# enkfsa.velocity <- lapply(enkfsa.velocity, function(s) lapply(s, as.matrix))

## Read bed and friction output from EnKFSA
bed_files <- lapply(1:length(sample_ind), function(s) paste0(enkfsa.output_dir[[s]], "/enkf_bed_sample", sample_ind[s], "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.bed <- lapply(bed_files, qread) #(file = paste0(enkfsa.output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
# enkfsa.bed_ini <- lapply(enkfsa.bed, function(s) s[[1]])
# enkfsa.bed_mid <- lapply(enkfsa.bed, function(s) s[[2]])
enkfsa.bed_fin <- lapply(enkfsa.bed, function(s) s[[3]])

fric_files <- lapply(1:length(sample_ind), function(s) paste0(enkfsa.output_dir[[s]], "/enkf_friction_sample", sample_ind[s], "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
enkfsa.fric <- lapply(fric_files, qread) #(file = paste0(enkfsa.output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date,  ".qs", sep = ""))
# enkfsa.fric_ini <- lapply(enkfsa.fric, function(s) s[[1]])
# enkfsa.fric_mid <- lapply(enkfsa.fric, function(s) s[[2]])
enkfsa.fric_fin <- lapply(enkfsa.fric, function(s) s[[3]])
# enkfsa.fric_ini <- lapply(enkfsa.fric_ini, exp)
# enkfsa.fric_mid <- lapply(enkfsa.fric_mid, exp)
enkfsa.fric_fin <- lapply(enkfsa.fric_fin, exp)

#################################
##      CRPS calculation       ##
#################################

crps_multivariate <- function(pred_samples, observed) {
  # Check dimensions
  n <- nrow(pred_samples)
  d <- ncol(pred_samples)
  
  # Ensure observed is a row vector
  if (length(observed) != d) {
    stop("Observed vector must have the same dimensions as the prediction samples.")
  }
  
  # Calculate the first term: Mean L2 norm between each sample vector and observed vector
  term1 <- mean(apply(pred_samples, 1, function(x) sqrt(sum((x - observed)^2))))
  
  # Calculate the second term: Mean L2 norm between all pairs of sample vectors
  term2 <- mean(apply(pred_samples, 1, function(x1) {
    mean(apply(pred_samples, 1, function(x2) sqrt(sum((x1 - x2)^2))))
  }))
  
  # Energy Score (Multivariate CRPS) is the difference between these two terms
  energy_score <- term1 - (term2 / 2)
  return(energy_score)
}


### To calculate RMSE for friction, need true gl for each sample
true_gl_pos <- true_gl[, years] * 1000 # Take GL position at last time point for each sample
true_gl_ind <- sapply(1:length(true_gl_pos), function(x) which(domain == true_gl_pos[x]))#floor(true_gl_pos / 800 * 2001) # Convert to grid index

if (recalculate_crps) {
    nsamples <- sapply(cnn.thickness_fin, function(x) dim(x)[2]) # minimum number of "particles", just in case there are some EnKF runs that fail
    state_samples <- sample(min(nsamples), 100) # use 1000 samples from the concatenated ensemble
    
    ## CRPS for the bed
    print("Calculating CRPS for CNN...")
    cnn.bed_crps <- c()
    cnn.fric_crps <- c()
    cnn.thickness_crps <- c()
    for (r in 1:length(s)) {
        
        cnn.bed_crps[r] <- crps_multivariate(t(cnn.bed_samples[[r]]), true_bed[r, ])
        cnn.fric_crps[r] <- crps_multivariate(t(cnn.fric_samples[[r]][1:true_gl_ind[r], ]), true_fric[r, 1:true_gl_ind[r]])
        cnn.thickness_crps[r] <- crps_multivariate(t(cnn.thickness_fin[[r]][, state_samples]), true_thicknesses_fin[r, ])
    }

    print("Calculating CRPS for Aug EnKF...")
    enkfsa.bed_crps <- c()
    enkfsa.fric_crps <- c()
    enkfsa.thickness_crps <- c()

    for (r in 1:length(s)) { # for each sample,
        enkfsa.bed_crps[r] <- crps_multivariate(t(enkfsa.bed_fin[[r]]), true_bed[r, ])
        enkfsa.fric_crps[r] <- crps_multivariate(t(enkfsa.fric_fin[[r]][1:true_gl_ind[r], ]), true_fric[r, 1:true_gl_ind[r]])
        enkfsa.thickness_crps[r] <- crps_multivariate(t(enkfsa.thickness_fin[[r]]), true_thicknesses_fin[r, ])
    }

    ## CRPS table
    cnn.bed_crps_mean <- mean(cnn.bed_crps)
    cnn.fric_crps_mean <- mean(cnn.fric_crps)
    cnn.thickness_crps_mean <- mean(cnn.thickness_crps)

    enkfsa.bed_crps_mean <- mean(enkfsa.bed_crps)
    enkfsa.fric_crps_mean <- mean(enkfsa.fric_crps)
    enkfsa.thickness_crps_mean <- mean(enkfsa.thickness_crps)

    crps_table <- data.frame(
        method = c("CNN", "EnKFSA"),
        thickness_fin = c(cnn.thickness_crps_mean, enkfsa.thickness_crps_mean),
        bed = c(cnn.bed_crps_mean, enkfsa.bed_crps_mean),
        friction = c(cnn.fric_crps_mean, enkfsa.fric_crps_mean)
    )
        write.csv(crps_table, file = paste0("./output/crps_table_", output_date, ".csv"))
        qsave(crps_table, file = paste0("./output/crps_table_", output_date, ".qs"))
    } else {
        crps_table <- qread(file = paste0("./output/crps_table_", output_date, ".qs"))
    }

print(crps_table)

#################################################################
## Univariate CRPS (per grid point) ##

crps_univ <- function(pred_samples, observed) { ## Univariate CRPS
  n <- length(pred_samples)
  
  # Calculate the first term: Average absolute error between samples and observed
  term1 <- mean(abs(pred_samples - observed))
  
  # Calculate the second term: Average absolute error between all pairs of samples
  term2 <- mean(abs(outer(pred_samples, pred_samples, "-")))
  
  # CRPS is the difference between these two terms
  crps <- term1 - (term2 / 2)
  return(crps)
}

print("Calculating CRPS per grid point...")
J <- dim(true_bed)[2]
cnn.bed_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)
cnn.fric_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)
cnn.thickness_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)
enkfsa.bed_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)
enkfsa.fric_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)
enkfsa.thickness_crps_per_gridpt_mat <- matrix(NA, nrow = length(s), ncol = J)

t1 <- system.time({
    for (r in 1:length(s)) {
    # cnn.bed_crps_per_gridpt <- c()
    cnn.fric_crps_per_gridpt <- c()
    # cnn.thickness_crps_per_gridpt <- c()
    # enkfsa.bed_crps_per_gridpt <- c()
    enkfsa.fric_crps_per_gridpt <- c()
    # enkfsa.thickness_crps_per_gridpt <- c()
    for (j in 1:J) {
        # cnn.bed_crps_per_gridpt[j] <- crps_univ(cnn.bed_samples[[r]][j, ], true_bed[r, j])
        cnn.fric_crps_per_gridpt[j] <- crps_univ(cnn.fric_samples[[r]][j, ], true_fric[r, j])
        # cnn.thickness_crps_per_gridpt[j] <- crps_univ(cnn.thickness_fin[[r]][j, ], true_thicknesses_fin[r, j])
        # enkfsa.bed_crps_per_gridpt[j] <- crps_univ(enkfsa.bed_fin[[r]][j, ], true_bed[r, j])
        enkfsa.fric_crps_per_gridpt[j] <- crps_univ(enkfsa.fric_fin[[r]][j, ], true_fric[r, j])
        # enkfsa.thickness_crps_per_gridpt[j] <- crps_univ(enkfsa.thickness_fin[[r]][j, ], true_thicknesses_fin[r, j])
    }
    # cnn.bed_crps_per_gridpt_mat[r, ] <- cnn.bed_crps_per_gridpt
    cnn.fric_crps_per_gridpt_mat[r, ] <- cnn.fric_crps_per_gridpt
    # cnn.thickness_crps_per_gridpt_mat[r, ] <- cnn.thickness_crps_per_gridpt
    # enkfsa.bed_crps_per_gridpt_mat[r, ] <- enkfsa.bed_crps_per_gridpt
    enkfsa.fric_crps_per_gridpt_mat[r, ] <- enkfsa.fric_crps_per_gridpt
    # enkfsa.thickness_crps_per_gridpt_mat[r, ] <- enkfsa.thickness_crps_per_gridpt

}

# cnn.mean_bed_crps_per_gridpt <- colMeans(cnn.bed_crps_per_gridpt_mat)
cnn.mean_fric_crps_per_gridpt <- colMeans(cnn.fric_crps_per_gridpt_mat)
# cnn.mean_thickness_crps_per_gridpt <- colMeans(cnn.thickness_crps_per_gridpt_mat)
# enkfsa.mean_bed_crps_per_gridpt <- colMeans(enkfsa.bed_crps_per_gridpt_mat)
enkfsa.mean_fric_crps_per_gridpt <- colMeans(enkfsa.fric_crps_per_gridpt_mat)
# enkfsa.mean_thickness_crps_per_gridpt <- colMeans(enkfsa.thickness_crps_per_gridpt_mat)

})

## Plot CRPS per grid point
if (save_plots) {
    print("Plotting CRPS per grid point...")
    # bed_crps_df <- data.frame(
    #     grid_point = 1:J,
    #     domain = domain/1000,
    #     cnn = cnn.mean_bed_crps_per_gridpt,
    #     enkfsa = enkfsa.mean_bed_crps_per_gridpt
    # )
    
    # bed_crps_plot <- ggplot(bed_crps_df) + 
    #     geom_line(aes(x = domain, y = enkfsa, color = "Aug EnKF"), lwd = 1) +
    #     geom_line(aes(x = domain, y = cnn, color = "CNN-EnKF"), lwd = 1) +
    #     labs(x = "Domain (km)", y = "CRPS", title = "Average CRPS per grid point (bed)") +
    #     theme_bw() +
    #     theme(text = element_text(size = 24)) +
    #     # theme(legend.position = "top") +
    #     scale_color_manual(values = c("turquoise", "red")) +
    #     theme(legend.title = element_blank())
    # png(paste0("./plots/bed_crps_per_gridpt_", output_date, ".png"), width = 1200, height = 500)
    # print(bed_crps_plot)
    
    # dev.off()
    
    min_GL <- min(true_gl_ind)
    fric_crps_df <- data.frame(
        grid_point = 1:min_GL,
        domain = domain[1:min_GL]/1000,
        cnn = cnn.mean_fric_crps_per_gridpt[1:min_GL],
        enkfsa = enkfsa.mean_fric_crps_per_gridpt[1:min_GL]
    )

    fric_crps_plot <- ggplot(fric_crps_df) + 
        geom_line(aes(x = domain, y = enkfsa, color = "Aug EnKF"), lwd = 1) +
        geom_line(aes(x = domain, y = cnn, color = "CNN-EnKF"), lwd = 1) +
        labs(x = "Domain (km)", y = "CRPS", title = "Average CRPS per grid point (friction)") +
        theme_bw() +
        theme(text = element_text(size = 24)) +
        # theme(legend.position = "top") +
        scale_color_manual(values = c("turquoise", "red")) +
        theme(legend.title = element_blank())

    png(paste0("./plots/fric_crps_per_gridpt_", output_date, ".png"), width = 1200, height = 500)
    print(fric_crps_plot)
    dev.off()

    # thickness_crps_df <- data.frame(
    #     grid_point = 1:J,
    #     domain = domain/1000,
    #     cnn = cnn.mean_thickness_crps_per_gridpt,
    #     enkfsa = enkfsa.mean_thickness_crps_per_gridpt
    # )

    # thickness_crps_plot <- ggplot(thickness_crps_df) + 
    #     geom_line(aes(x = domain, y = enkfsa, color = "Aug EnKF"), lwd = 1) +
    #     geom_line(aes(x = domain, y = cnn, color = "CNN-EnKF"), lwd = 1) +
    #     labs(x = "Domain (km)", y = "CRPS", title = "Average CRPS per grid point (thickness)") +
    #     theme_bw() +
    #     theme(text = element_text(size = 24)) +
    #     # theme(legend.position = "top") +
    #     scale_color_manual(values = c("turquoise", "red")) +
    #     theme(legend.title = element_blank())

    # png(paste0("./plots/thickness_crps_per_gridpt_", output_date, ".png"), width = 1200, height = 500)
    # print(thickness_crps_plot)
    # dev.off()
}
