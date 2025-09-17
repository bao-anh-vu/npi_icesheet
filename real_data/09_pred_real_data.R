## Predict on real data

## Post-process output from CNN
setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())
source("./source/create_model.R")
source("./source/posterior_loss.R")
source("./source/sample_from_posterior.R")

# library(parallel)
library(mvtnorm) 
library(Matrix)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(qs)
library(abind)
library(dplyr)

## Flags
save_pred <- T
save_plots <- T
log_transform <- T
test_on_train <- F
use_missing_pattern <- T

## Read data
data_date <- "20241111" #"20241103"
sets <- 1:20 #6:20
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

if (use_missing_pattern) {
    data_dir <- paste0("./data/training_data/", setsf, "/missing/")
    output_dir <- paste0("./output/cnn/", setsf, "/missing/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/missing/")
} else {
    data_dir <- paste0("./data/training_data/", setsf, "/nonmissing/")
    output_dir <- paste0("./output/cnn/", setsf, "/nonmissing/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/nonmissing/")
}

ssa_steady <- qread(file = paste0("data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain

test_data <- qread(file = paste0(data_dir, "/test_data_", data_date, ".qs"))

## Also load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
# velocity_data <- qread(file = "./data/velocity/all_velocity_arr.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

## Replace NA values in real_data with 0
surf_elev_data[is.na(surf_elev_data)] <- 0
velocity_data[is.na(velocity_data)] <- 0

input_mean <- test_data$input_mean
input_sd <- test_data$input_sd
mean_surf_elev <- input_mean[1] #mean(surf_elev_data, na.rm = T)
sd_surf_elev <- input_sd[1] #sd(surf_elev_data, na.rm = T)
surf_elev_data_std <- (surf_elev_data - mean_surf_elev) / sd_surf_elev

mean_velocity <- input_mean[2] #mean(velocity_data, na.rm = T)
sd_velocity <- input_sd[2] #sd(velocity_data, na.rm = T)
velocity_data_std <- (velocity_data - mean_velocity) / sd_velocity

real_data <- abind(surf_elev_data_std, velocity_data_std, along = 3)
real_data <- array(real_data, dim = c(1L, dim(real_data)))

## Have to standardise the input
test_input <- test_data$input
test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)

sim <- sample(1:dim(test_data$input)[1], 1)
png("./plots/cnn/test_vs_real.png", width = 1000, height = 800)

par(mfrow = c(2,1))
matplot(test_data$input[sim,,,1], col = "grey", type = "l", main = paste0("Simulation ", sim))
matlines(real_data[1,,,1], col = "salmon")

matplot(test_data$input[sim,,,2], col = "grey", type = "l", main = paste0("Simulation ", sim))
matlines(real_data[1,,,2], col = "salmon")
dev.off()

input_dim <- dim(test_data$input)[2:4]
n_basis_funs <- dim(test_data$fric_coefs)[2]
n_gl <- dim(test_data$grounding_line)[2]
n_mean_elements <- n_basis_funs * 2 + n_gl
n_chol_elements <- (n_basis_funs * 2 + n_gl) + (n_basis_funs - 1) * 2 + (n_gl - 1) # diagonal + lower-diag elements
output_dim <- n_mean_elements + n_chol_elements  # THIS NEEDS TO CHANGE TO THE TOTAL NUMBER OF BASIS FUNCTIONS + COVARIANCE PARAMETERS


## Load the model
model <- create_model_posterior(input_dim = input_dim, 
                                output_dim = output_dim,
                                n_basis_funs = n_basis_funs,
                                n_gl = n_gl)

### Display the model's architecture
summary(model)

### Plot the loss
history <- qread(file = paste0(output_dir, "history_", data_date, ".qs"))
# history %>% plot() #+
# coord_cartesian(xlim = c(1, epochs))

# if (use_missing_pattern) {
#     plot_dir <- paste0("./plots/cnn/", setsf, "/missing")
# } else {
#     plot_dir <- paste0("./plots/cnn/", setsf, "/nonmissing")

# }

# if (save_plots) {

#     if (!dir.exists(plot_dir)) {
#         dir.create(plot_dir)
#     }

#     png(paste0(plot_dir, "/loss.png"), width = 1000, height = 500)
#     plot(history$metrics$loss[2:100], type = "l")
#     lines(history$metrics$val_loss[2:100], col = "red")
#     legend("topright", legend = c("Training", "Validation"), col = c("black", "red"), lty = 1, cex = 0.8)
#     dev.off()
# }

###################################
##        Mean prediction        ##   
###################################

## Reload model from checkpoint
checkpoint_path <- paste0(output_dir, "/checkpoints/cp-{epoch:04d}.ckpt")
checkpoint_dir <- fs::path_dir(checkpoint_path)

checkpt <- which.min(history$metrics$val_loss) 
if (!is.null(checkpt)) {
    ## Load the model from checkpt
    cp_restart <- paste0(output_dir, "/checkpoints/cp-", formatC(checkpt, width = 4, flag = "0"), ".ckpt")
    # latest <- tf$train$latest_checkpoint(checkpoint_dir)
    load_model_weights_tf(model, cp_restart)
} else { # use latest checkpt
    latest <- tf$train$latest_checkpoint(checkpoint_dir)
    load_model_weights_tf(model, latest)
}

## Predict
if (test_on_train) {
    train_subset <- train_data$input[1:100,,,]
    pred_output <- model %>% predict(train_subset)
} else {

    pred_time <- system.time({
        # pred_output <- model %>% predict(test_input)
        pred_output <- model %>% predict(real_data)

    })
    cat("Prediction time: ", pred_time[3], "\n")
}


## Plot the predicted friction coefs generated by these basis functions

fric_basis_mat <- test_data$fric_basis_mat
bed_basis_mat <- test_data$bed_basis_mat
n_fric_basis <- ncol(fric_basis_mat)
n_bed_basis <- ncol(bed_basis_mat)

pred_mean <- pred_output[, 1:n_mean_elements]

## Un-standardise output
pred_fric_coefs <- pred_mean[1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs <- pred_mean[(n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl <- pred_mean[(n_fric_basis+n_bed_basis+1):length(pred_mean)] * test_data$sd_gl + test_data$mean_gl

# if (test_on_train) {
#     train_fric_coefs <- train_data$fric_coefs[1:100, ] * train_data$sd_fric_coefs + train_data$mean_fric_coefs
# }
test_fric_coefs <- test_output[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
test_bed_coefs <- test_output[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
test_gl <- test_output[, (n_fric_basis+n_bed_basis+1):ncol(test_output)] * test_data$sd_gl + test_data$mean_gl

png(paste0(plot_dir, "/pred_fric_real.png"), width = 2000, height = 1200)
plot(pred_fric_coefs, type = "l")
# lines(test_fric_coefs[1,], col = "red")
dev.off()

png(paste0(plot_dir, "/pred_bed_real.png"), width = 2000, height = 1200)
plot(pred_bed_coefs, type = "l")
# lines(test_bed_coefs[1,], col = "red")
dev.off()

## Compute predicted vs original friction 
if (log_transform) {
    pred_fric <- exp(fric_basis_mat %*% pred_fric_coefs)
    test_fric <- exp(fric_basis_mat %*% t(test_fric_coefs))
} else {
    pred_fric <- fric_basis_mat %*% pred_fric_coefs
    test_fric <- fric_basis_mat %*% t(test_fric_coefs)
}

## Compute predicted vs original bed fluctuations (demeaned bed)
pred_bed_demean <- bed_basis_mat %*% pred_bed_coefs
test_bed_demean <- bed_basis_mat %*% t(test_bed_coefs)

# if (test_on_train) {
#     par(mfrow = c(1,1))
#     plot(train_fric_coefs[50,], type = "l")
#     lines(pred_fric_coefs[50,], col = "red")
# } else {
#     par(mfrow = c(1,1))
#     plot(test_fric_coefs[1,], type = "l")
#     lines(pred_fric_coefs[1,], col = "red")
# }

### Reconstruct bed elevation by adding bed trend to the predicted oscillations
# bed_obs <- readRDS(file = paste0("./data/training_data/bed_obs_", data_date, ".rds"))

bed_basis <- qread(file = paste0("./data/training_data/bed_basis_01_", data_date, ".qs"))

bed_mean <- bed_basis$mean
# bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = ncol(pred_bed_demean))
pred_bed <- pred_bed_demean + bed_mean
# test_bed <- test_bed_demean + bed_mean_mat

## Covariance matrix
pred_chol <- pred_output[, (n_mean_elements+1):ncol(pred_output)]

## Construct Cholesky factor of the precision
Lb_elems <- n_basis_funs + (n_basis_funs - 1)
Lc_elems <- n_basis_funs + (n_basis_funs - 1)
Lg_elems <- n_gl + (n_gl - 1)

Lmats <- list()
# for (s in 1:nrow(pred_chol)) {
    Lb <- construct_L_matrix(pred_chol[1:Lb_elems], n_basis_funs)
    Lc <- construct_L_matrix(pred_chol[(Lb_elems+1):(Lb_elems+Lc_elems)], n_basis_funs)
    Lg <- construct_L_matrix(pred_chol[(Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)
    Lmat <- bdiag(Lb, Lc, Lg)
    # Lmats[[s]] <- bdiag(Lb, Lc, Lg)
# }

# ## Need to sample from the posterior distribution of the coefs
# ## then transform them to actual friction, bed, gl

S <- 1000 ## number of posterior samples

## Sample basis function coefficients from posterior
pred_samples_ls <- sample_from_posterior(n = S, mean = pred_mean, prec_chol = Lmat)    

## Transform basis coefficients into actual friction coefficients
sd_fric_coefs <- test_data$sd_fric_coefs
mean_fric_coefs <- test_data$mean_fric_coefs

fric_coefs <- pred_samples_ls[1:n_fric_basis, ]
fric_coefs_ustd <- fric_coefs * sd_fric_coefs + mean_fric_coefs
fric_samples <- exp(fric_basis_mat %*% fric_coefs_ustd)

## Transform basis coefficients into bed elevations
bed_coefs <- pred_samples_ls[(n_fric_basis+1):(n_fric_basis+n_bed_basis), ]
bed_coefs_ustd <- bed_coefs * test_data$sd_bed_coefs + test_data$mean_bed_coefs
bed_samples <- bed_basis_mat %*% bed_coefs_ustd + bed_mean

## Grounding line is a direct output from the CNN
gl_samples_ls <- pred_samples_ls[(n_fric_basis+n_bed_basis+1):n_mean_elements, ]

## Compute quantiles
fric_q <- apply(fric_samples, 1, quantile, probs = c(0.025, 0.975))
bed_q <- apply(bed_samples, 1, quantile, probs = c(0.025, 0.975))
gl_q <- apply(gl_samples_ls, 1, quantile, probs = c(0.025, 0.975))

fric_lq <- fric_q[1,]
fric_uq <- fric_q[2,]
bed_lq <- bed_q[1,]
bed_uq <- bed_q[2,]
gl_lq <- gl_q[1,]
gl_uq <- gl_q[2,]

## Plot the results
png(paste0(plot_dir, "/fric_samples.png"), width = 1000, height = 500)
plot(domain/1e3, fric_lq, type = "l", col = "grey", lwd = 2)
lines(domain/1e3, fric_uq, col = "grey", lwd = 2)
lines(domain/1e3, pred_fric, col = "red", lwd = 2)
dev.off()

## Validate against the rest of the bed obs
bed_obs_df <- qread(file = "./data/bed_obs_df.qs")
bed_obs_train <- bed_obs_df %>% filter(chosen == 1)
bed_obs_val <- bed_obs_df %>% filter(chosen == 0)

png(paste0(plot_dir, "/bed_samples.png"), width = 1000, height = 500)
plot(domain/1e3, bed_lq, type = "l", col = "grey", lwd = 2, xlab = "Grid point", ylab = "Elevation (m)")
lines(domain/1e3, bed_uq, col = "grey", lwd = 2)
lines(domain/1e3, pred_bed, col = "red", lwd = 2)
points(domain/1e3, bed_obs_train$ind, bed_obs_train$bed_elev, col = "black", pch = 20)
points(domain/1e3, bed_obs_val$ind, bed_obs_val$bed_elev, col = "cyan")
dev.off()

## Save predictions
if (save_pred) {
    qsave(pred_fric, file = paste0(output_dir, "/pred_fric_real_", data_date, ".qs"))
    qsave(pred_bed, file = paste0(output_dir, "/pred_bed_real_", data_date, ".qs"))
    qsave(fric_samples, file = paste0(output_dir, "/fric_samples_real_", data_date, ".qs"))
    qsave(bed_samples, file = paste0(output_dir, "/bed_samples_real_", data_date, ".qs"))
    qsave(fric_q, file = paste0(output_dir, "/fric_quantiles_real_", data_date, ".qs"))
    qsave(bed_q, file = paste0(output_dir, "/bed_quantiles_real_", data_date, ".qs"))
    qsave(gl_q, file = paste0(output_dir, "/gl_quantiles_real_", data_date, ".qs"))
}