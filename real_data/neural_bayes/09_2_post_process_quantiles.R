## Post-process output from CNN
setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())
source("./source/create_model.R")
source("./source/custom_loss_function.R")

library(mvtnorm) 
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(qs)
library(abind)
library(dplyr)

## Flags
# output_var <- "all" # "friction" # "grounding_line" # "bed_elevation
# sim_beds <- T
save_pred <- T
save_plots <- T
log_transform <- T
use_missing_pattern <- T
test_on_train <- F

## Read data
data_date <- "20241111"
sets <- 1:20 #6:20
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

# if (sim_beds) {
if (use_missing_pattern) {
    data_dir <- paste0("./data/training_data", "/", setsf, "/missing/")
} else {
    data_dir <- paste0("./data/training_data", "/", setsf, "/nonmissing/")
}

test_data <- qread(file = paste0(data_dir, "test_data_", data_date, ".qs"))

if (test_on_train) {
    test_data <- qread(file = paste0(data_dir, "train_data_", data_date, ".qs"))
}
test_input <- test_data$input
test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)

## Also load real data and standardise it to make it suitable for the CNN
# surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
surf_elev_data <- qread(file = paste0("./data/surface_elev/adj_se_mat_", data_date, ".qs"))
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


## Flowline data
ssa_steady <- qread(file = paste("./data/training_data/steady_state/steady_state_", data_date, ".qs", sep = ""))
domain <- ssa_steady$domain

## Bed observations
bed_obs_df <- qread(file = "./data/bed_obs_df.qs")
bed_obs_validate <- bed_obs_df %>% filter(chosen == 0 & !is.na(bed_elev) & !is.na(bed_sd))
bed_obs_chosen <- bed_obs_df %>% filter(chosen == 1)

##############################
##      Load the model      ##
##############################

if (test_on_train) {
    test_input <- test_data$input[1:100,,,]
    test_output <- test_output[1:100, ]
}

input_dim <- dim(test_input)[2:4]
output_dim <- ncol(test_output)
model <- create_model(input_dim = input_dim, output_dim = output_dim)

# Display the model's architecture
summary(model)

## Basis function info
friction_basis <- qread(file = paste0("data/training_data/friction_basis_01_", data_date, ".qs"))
fric_basis_mat <- friction_basis$basis_mat
n_fric_basis <- ncol(fric_basis_mat)

bed_basis <- qread(file = paste0("data/training_data/bed_basis_01_", data_date, ".qs"))
bed_basis_mat <- bed_basis$basis_mat
n_bed_basis <- ncol(bed_basis_mat)



# plot(history$metrics$loss, type = "l")
# lines(history$metrics$val_loss, col = "red")
# legend("topright", legend = c("Training", "Validation"), col = c("black", "red"), lty = 1, cex = 0.8)

# plot_dir <- paste0("./plots/", output_var, "/", setsf)

# if (save_plots) {

#     if (!dir.exists(plot_dir)) {
#         dir.create(plot_dir)
#     }



###################################
##        Mean prediction        ##   
###################################

if (use_missing_pattern) {
    output_dir_mean <- paste0("./output/neural_bayes/", setsf, "/missing/mean/")
    # plot_dir_mean <- paste0("./plots/neural_bayes/", setsf, "/missing/")
} else {
    output_dir_mean <- paste0("./output/neural_bayes/", setsf, "/nonmissing/mean/")
    # plot_dir_mean <- paste0("./plots/neural_bayes/", setsf, "/nonmissing/")
}

plot_dir <- paste0("./plots/neural_bayes/", setsf, "/")

history_mean <- qread(file = paste0(output_dir_mean, "history_", data_date, ".qs"))

# png(paste0(plot_dir_mean, "loss_mean_", data_date, ".png"), width = 1000, height = 500)
# plot(history_mean$metrics$loss[2:50], type = "l")
# lines(history_mean$metrics$val_loss[2:50], col = "red")
# legend("topright", legend = c("Training", "Validation"), col = c("black", "red"), lty = 1, cex = 0.8)
# dev.off()

## Reload model from checkpoint
checkpoint_path <- paste0(output_dir_mean, "checkpoints/cp-{epoch:04d}.ckpt")
checkpoint_dir <- fs::path_dir(checkpoint_path)

checkpt <- which.min(history_mean$metrics$val_loss) # c(69)
# for (checkpt in checkpts) {
if (!is.null(checkpt)) {
    ## Load the model from checkpt
    cp_restart <- paste0(output_dir_mean, "checkpoints/cp-", formatC(checkpt, width = 4, flag = "0"), ".ckpt")
    # latest <- tf$train$latest_checkpoint(checkpoint_dir)
    load_model_weights_tf(model, cp_restart)
} else {
    latest <- tf$train$latest_checkpoint(checkpoint_dir)
    load_model_weights_tf(model, latest)
}

## Predict
pred_output <- model %>% predict(test_input)
pred_output_real <- model %>% predict(real_data)

## Un-standardise output
pred_fric_coefs <- pred_output[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs <- pred_output[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl <- pred_output[, (n_fric_basis+n_bed_basis+1):ncol(pred_output)] * test_data$sd_gl + test_data$mean_gl

test_fric_coefs <- test_output[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
test_bed_coefs <- test_output[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
test_gl <- test_output[, (n_fric_basis+n_bed_basis+1):ncol(test_output)] * test_data$sd_gl + test_data$mean_gl

## Do the same for output on real data
pred_fric_coefs_real <- pred_output_real[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs_real <- pred_output_real[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl_real <- pred_output_real[, (n_fric_basis+n_bed_basis+1):ncol(pred_output_real)] * test_data$sd_gl + test_data$mean_gl

## Compute predicted vs original friction -- this needs more work!!
if (log_transform) {
    pred_fric <- exp(fric_basis_mat %*% t(pred_fric_coefs))
    test_fric <- exp(fric_basis_mat %*% t(test_fric_coefs))
    pred_fric_real <- exp(fric_basis_mat %*% pred_fric_coefs_real)
} else {
    pred_fric <- fric_basis_mat %*% t(pred_fric_coefs)
    test_fric <- fric_basis_mat %*% t(test_fric_coefs)
    pred_fric_real <- fric_basis_mat %*% pred_fric_coefs_real
}

## Compute predicted vs original bed fluctuations (demeaned bed)
pred_bed_demean <- bed_basis_mat %*% t(pred_bed_coefs)
test_bed_demean <- bed_basis_mat %*% t(test_bed_coefs)
pred_bed_demean_real <- bed_basis_mat %*% pred_bed_coefs_real

# if (test_on_train) {
#     par(mfrow = c(1,1))
#     plot(train_fric_coefs[50,], type = "l")
#     lines(pred_fric_coefs[50,], col = "red")
# } else {
#     par(mfrow = c(1,1))
#     plot(test_fric_coefs[1,], type = "l")
#     lines(pred_fric_coefs[1,], col = "red")
# }

############################
##         Median         ##
############################

output_dir_med <- paste0("./output/neural_bayes/", setsf, "/missing/q05/")
# plot_dir_quantiles <- paste0("./plots/cnn/neural_bayes/", setsf, "/")    
model_med <- create_model_quantile(input_dim = input_dim, output_dim = output_dim, quantile = 0.5)

history_med <- qread(file = paste0(output_dir_med, "history_", data_date, ".qs"))
# history %>% plot() #+
# coord_cartesian(xlim = c(1, epochs))
png(paste0(plot_dir, "/loss_median_", data_date, ".png"), width = 1000, height = 500)
# par(mfrow = c(2,1))
plot(history_med$metrics$loss, type = "l", 
        xlab = "Iteration", ylab = "Loss",
        main = paste0("Loss when training the 50th quantile"))
lines(history_med$metrics$val_loss, col = "red")
legend("topright", legend = c("Training", "Validation"), 
        col = c("black", "red"), lty = 1, cex = 0.8)
dev.off()

## Reload model for quantile from checkpoint
checkpoint_path_med <- paste0(output_dir_med, "/checkpoints/cp-{epoch:04d}.ckpt")
checkpoint_dir_med <- fs::path_dir(checkpoint_path_med)

checkpt_med <- which.min(history_med$metrics$val_loss) # c(69)
# for (checkpt in checkpts) {
if (!is.null(checkpt_med)) {
    ## Load the model from checkpt
    cp_restart_med <- paste0(output_dir_med, "/checkpoints/cp-", formatC(checkpt_med, width = 4, flag = "0"), ".ckpt")
    # latest <- tf$train$latest_checkpoint(checkpoint_dir)
    load_model_weights_tf(model_med, cp_restart_med)
} else {
    latest_cp_med <- tf$train$latest_checkpoint(checkpoint_dir_med)
    load_model_weights_tf(model_med, latest_cp_med)
}

## Predict
pred_output_med <- model_med %>% predict(test_input)

## Predict on real data
pred_output_med_real <- model_med %>% predict(real_data)

## Un-standardise output
pred_fric_coefs_med <- pred_output_med[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs_med <- pred_output_med[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl_med <- pred_output_med[, (n_fric_basis+n_bed_basis+1):ncol(pred_output_med)] * test_data$sd_gl + test_data$mean_gl

test_fric_coefs <- test_output[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
test_bed_coefs <- test_output[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
test_gl <- test_output[, (n_fric_basis+n_bed_basis+1):ncol(test_output)] * test_data$sd_gl + test_data$mean_gl

## Do the same for output on real data
pred_fric_coefs_med_real <- pred_output_med_real[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs_med_real <- pred_output_med_real[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl_med_real <- pred_output_med_real[, (n_fric_basis+n_bed_basis+1):ncol(pred_output_med_real)] * test_data$sd_gl + test_data$mean_gl

# ######################################
# ##        Quantile prediction       ##
# ######################################
quantiles <- c(0.025, 0.975)
pred_quantiles <- list()
pred_quantiles_real <- list()

for (q in 1:length(quantiles)) {
    quantile <- quantiles[q]
    q_string <- gsub("\\.", "", as.character(quantile))
    output_dir_quantile <- paste0("./output/neural_bayes/", setsf, "/missing/q", q_string, "/")
    
    model_quantile <- create_model_quantile(input_dim = input_dim, output_dim = output_dim, quantile = quantile)

    history_quantile <- qread(file = paste0(output_dir_quantile, "history_", data_date, ".qs"))
    # history %>% plot() #+
    # coord_cartesian(xlim = c(1, epochs))
    png(paste0(plot_dir, "/loss_quantile_", q_string, "_", data_date, ".png"), width = 1000, height = 500)
    par(mfrow = c(1,1))
    plot(history_quantile$metrics$loss, type = "l", 
            xlab = "Iteration", ylab = "Loss",
            main = paste0("Loss when training the", quantile*100, "th quantile"))
    lines(history_quantile$metrics$val_loss[2:100], col = "red")
    legend("topright", legend = c("Training", "Validation"), 
            col = c("black", "red"), lty = 1, cex = 0.8)
    dev.off()
    
    ## Reload model for quantile from checkpoint
    checkpoint_path_quantile <- paste0(output_dir_quantile, "checkpoints/cp-{epoch:04d}.ckpt")
    checkpoint_dir_quantile <- fs::path_dir(checkpoint_path_quantile)

    checkpt_quantile <- which.min(history_quantile$metrics$val_loss) # c(69)
    # for (checkpt in checkpts) {
    if (!is.null(checkpt_quantile)) {
        ## Load the model from checkpt
        cp_restart_quantile <- paste0(output_dir_quantile, "checkpoints/cp-", formatC(checkpt_quantile, width = 4, flag = "0"), ".ckpt")
        # latest <- tf$train$latest_checkpoint(checkpoint_dir)
        load_model_weights_tf(model_quantile, cp_restart_quantile)
    } else {
        latest_cp_quantile <- tf$train$latest_checkpoint(checkpoint_dir_quantile)
        load_model_weights_tf(model_quantile, latest_cp_quantile)
    }

    ## Predict
    pred_output_quantile <- model_quantile %>% predict(test_input)

    ## Predict on real data
    pred_output_quantile_real <- model_quantile %>% predict(real_data)

    ## Un-standardise output
    pred_fric_coefs_quantile <- pred_output_quantile[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
    pred_bed_coefs_quantile <- pred_output_quantile[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
    pred_gl_quantile <- pred_output_quantile[, (n_fric_basis+n_bed_basis+1):ncol(pred_output_quantile)] * test_data$sd_gl + test_data$mean_gl

    ## Do the same for output on real data
    pred_fric_coefs_quantile_real <- pred_output_quantile_real[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
    pred_bed_coefs_quantile_real <- pred_output_quantile_real[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
    pred_gl_quantile_real <- pred_output_quantile_real[, (n_fric_basis+n_bed_basis+1):ncol(pred_output_quantile_real)] * test_data$sd_gl + test_data$mean_gl

    # pred_quantiles[[q]] <- list(fric = pred_fric_quantile, bed = pred_bed_quantile, gl = pred_gl_quantile)
    pred_quantiles[[q]] <- list(fric = pred_fric_coefs_quantile, 
                                bed = pred_bed_coefs_quantile, 
                                gl = pred_gl_quantile)

    pred_quantiles_real[[q]] <- list(fric = pred_fric_coefs_quantile_real, 
                                     bed = pred_bed_coefs_quantile_real, 
                                     gl = pred_gl_quantile_real)
    
}

if (test_on_train) {
    plot_tag <- "_train"
} else {
    plot_tag <- "_test"
}

s <- 1
png(paste0("./plots/neural_bayes/", setsf, "/pred_coefs", plot_tag, ".png"), width = 1000, height = 1000)
par(mfrow = c(2,1))
plot(test_fric_coefs[s,], lwd = 1.5, #ylim = c(-5, 5),
    type = "l", main = "Friction basis coefficients")
lines(pred_fric_coefs[s,], col = "red")
# lines(pred_fric_coefs_med[s,], col = "red")
lines(pred_quantiles[[1]]$fric[s, ], col = "salmon")
lines(pred_quantiles[[2]]$fric[s, ], col = "salmon")

plot(test_bed_coefs[s,], lwd = 1.5, type = "l", 
    main = "Bed elevation basis coefficients")
lines(pred_bed_coefs[s,], col = "red")
# lines(pred_bed_coefs_med[s,], col = "red")
lines(pred_quantiles[[1]]$bed[s, ], col = "salmon")
lines(pred_quantiles[[2]]$bed[s, ], col = "salmon")
dev.off()

## Plot predicted coefficients on real data too
png(paste0("./plots/neural_bayes/", setsf, "/pred_coefs_real.png"), width = 1000, height = 1000)
par(mfrow = c(2,1))
plot(pred_fric_coefs_real, lwd = 1.5, col = "red", #ylim = c(-5, 5),
    type = "l", main = "Friction basis coefficients (real data)")
# lines(pred_fric_coefs_med_real, col = "blue")
lines(pred_quantiles_real[[1]]$fric, col = "salmon")
lines(pred_quantiles_real[[2]]$fric, col = "salmon")

plot(pred_bed_coefs_real, lwd = 1.5, col = "red", type = "l", 
    main = "Bed elevation basis coefficients (real data)")
# lines(pred_bed_coefs_med_real, col = "blue")
lines(pred_quantiles_real[[1]]$bed, col = "salmon")
lines(pred_quantiles_real[[2]]$bed, col = "salmon")
dev.off()

# ########### end quantile ################

## Calculate quantiles of the bed and friction based on 
## the quantiles of the basis function coefficients

### Friction
z_val <- qnorm(0.975)
fric_sd1 <- (pred_fric_coefs - pred_quantiles[[1]]$fric) / z_val # distance from mean to lower quantile
fric_sd2 <- (pred_quantiles[[2]]$fric - pred_fric_coefs) / z_val # distance from mean to upper quantile
fric_avg_sd <- (fric_sd1 + fric_sd2)/2 # this is the "point-wise" standard deviation 

neg_ind <- which(fric_avg_sd < 0)
fric_avg_sd[neg_ind] <- 1e-02
## Now sample lots of basis coefficients from this distribution
## Maybe for each test output sample, we can sample 1000 basis coefficients
n_fric_basis <- ncol(fric_basis_mat)
fric_coef_samples <- list()
fric_lq <- list()
fric_uq <- list()

S <- 1000
# r <- 213
for (r in 1:nrow(test_output)) { ## need to parallelise this so that it's faster
    
    fric_coef_samples_list <- lapply(1:S, function(i) rnorm(n_fric_basis, 
                                                            mean = pred_fric_coefs[r, ], 
                                                            sd = fric_avg_sd[r, ]))

    fric_coef_samples <- t(matrix(unlist(fric_coef_samples_list), nrow = n_fric_basis, ncol = S))
    
    ## Exponentiate here to get friction on original scale
    if (log_transform) {
        fric_samples <- exp(fric_basis_mat %*% t(fric_coef_samples))
    } else {
        fric_samples <- fric_basis_mat %*% t(fric_coef_samples)
    }

    ## Now take the pointwise quantiles
    fric_lq[[r]] <- apply(fric_samples, 1, quantile, probs = 0.025)
    fric_uq[[r]] <- apply(fric_samples, 1, quantile, probs = 0.975)
    
}

## Do the same for real data
fric_sd1_real <- (pred_fric_coefs_real - pred_quantiles_real[[1]]$fric) / z_val # distance from mean to lower quantile
fric_sd2_real <- (pred_quantiles_real[[2]]$fric - pred_fric_coefs_real) / z_val # distance from mean to upper quantile
fric_avg_sd_real <- (fric_sd1 + fric_sd2)/2 # this is the "point-wise" standard deviation 

fric_coef_samples_list_real <- lapply(1:S, function(i) rnorm(n_fric_basis, 
                                                        mean = pred_fric_coefs_real, 
                                                        sd = fric_avg_sd_real))
fric_coef_samples_real <- t(matrix(unlist(fric_coef_samples_list_real), nrow = n_fric_basis, ncol = S))
    
## Exponentiate here to get friction on original scale
if (log_transform) {
    fric_samples_real <- exp(fric_basis_mat %*% t(fric_coef_samples_real))
} else {
    fric_samples_real <- fric_basis_mat %*% t(fric_coef_samples_real)
}

## Now take the pointwise quantiles
fric_lq_real <- apply(fric_samples_real, 1, quantile, probs = 0.025)
fric_uq_real <- apply(fric_samples_real, 1, quantile, probs = 0.975)

### Bed elevation

### Add bed trend to the predicted oscillations
# bed_obs <- qread(file = paste0("./training_data/bed_obs_", data_date, ".qs"))
# df <- data.frame(obs_locations = domain[bed_obs$locations], bed_elev = bed_obs$obs)
# bed.fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.25, 
#                     control = loess.control(surface = "direct")) 
# bed_mean <- predict(bed.fit, newdata = data.frame(obs_locations = domain))
bed_mean <- bed_basis$mean
bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = ncol(pred_bed_demean))

pred_bed <- pred_bed_demean + bed_mean_mat
test_bed <- test_bed_demean + bed_mean_mat
pred_bed_real <- pred_bed_demean_real + bed_mean

#### Quantiles for bed oscillations
z_val <- qnorm(0.975)
bed_sd1 <- (pred_bed_coefs - pred_quantiles[[1]]$bed) / z_val # distance from mean to lower quantile
bed_sd2 <- (pred_quantiles[[2]]$bed - pred_bed_coefs) / z_val # distance from mean to upper quantile
bed_avg_sd <- (bed_sd1 + bed_sd2)/2 # this is the "point-wise" standard deviation 

## Now sample lots of basis coefficients from this distribution
## Maybe for each test output sample, we can sample 1000 basis coefficients
n_bed_basis <- ncol(bed_basis_mat)
bed_coef_samples <- list()
bed_lq <- list()
bed_uq <- list()

bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = S)
for (r in 1:nrow(test_output)) {

    cat("r =", r, "\n")
    bed_coef_samples_list <- lapply(1:S, function(i) try(rnorm(n_bed_basis, 
                                                                mean = pred_bed_coefs[r, ], 
                                                                sd = bed_avg_sd[r, ]))
                                    )


    

    bed_coef_samples <- t(matrix(unlist(bed_coef_samples_list), nrow = n_bed_basis, ncol = S))
    bed_samples_demean <- bed_basis_mat %*% t(bed_coef_samples)
    bed_samples <- bed_samples_demean + bed_mean_mat

    ## Now take the pointwise quantiles
    bed_lq[[r]] <- apply(bed_samples, 1, quantile, probs = 0.025)
    bed_uq[[r]] <- apply(bed_samples, 1, quantile, probs = 0.975)
    
}

## Do the same for real data
bed_sd1_real <- (pred_bed_coefs_real - pred_quantiles_real[[1]]$bed) / z_val # distance from mean to lower quantile
bed_sd2_real <- (pred_quantiles_real[[2]]$bed - pred_bed_coefs_real) / z_val # distance from mean to upper quantile
bed_avg_sd_real <- (bed_sd1 + bed_sd2)/2 # this is the "point-wise" standard deviation

bed_coef_samples_list_real <- lapply(1:S, function(i) try(rnorm(n_bed_basis, 
                                                            mean = pred_bed_coefs_real, 
                                                            sd = bed_avg_sd_real))
                                )
bed_coef_samples_real <- t(matrix(unlist(bed_coef_samples_list_real), nrow = n_bed_basis, ncol = S))
bed_samples_demean_real <- bed_basis_mat %*% t(bed_coef_samples_real)
bed_samples_real <- bed_samples_demean_real + bed_mean_mat

bed_lq_real <- apply(bed_samples_real, 1, quantile, probs = 0.025)
bed_uq_real <- apply(bed_samples_real, 1, quantile, probs = 0.975)

### Grounding line
pred_gl_lq <- pred_quantiles[[1]]$gl
pred_gl_uq <- pred_quantiles[[2]]$gl

#################################
##      Save predictions       ##
#################################

if (save_pred) {
    qsave(pred_fric, file = paste0(output_dir_mean, "/pred_fric_", data_date, ".qs"))
    qsave(pred_bed, file = paste0(output_dir_mean, "/pred_bed_", data_date, ".qs"))
    qsave(pred_gl, file = paste0(output_dir_mean, "/pred_gl_", data_date, ".qs"))
    qsave(pred_fric_real, file = paste0(output_dir_mean, "/pred_fric_real_", data_date, ".qs"))
    qsave(pred_bed_real, file = paste0(output_dir_mean, "/pred_bed_real_", data_date, ".qs"))
    qsave(pred_gl_real, file = paste0(output_dir_mean, "/pred_gl_real_", data_date, ".qs"))
}


######################################
##      Comparison with truth       ##
######################################

## True parameter values for comparison
true_fric <- test_data$true_fric
true_bed <- test_data$true_bed
true_gl <- test_data$true_gl #* test_data$sd_gl + test_data$mean_gl

## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1/3)
# true_fric <- true_fric/fric_scale

## Choose random samples from test set for plotting
samples <- sample(1:nrow(test_output), 4)

if (save_plots) {
    
    ## Friction plots
    png(paste0(plot_dir, "/pred_vs_true_fric.png"), width = 2000, height = 1200)

    par(mfrow = c(length(samples) / 2, 2))
    for (i in 1:length(samples)) {
        par(mar = c(6, 8, 4, 2))
        # gl_i <- test_data$grounding_line[i] / domain[length(domain)] * length(domain)
        plot_domain <- 1:length(domain) # ceiling(gl)

        plot(domain[plot_domain]/1000, test_fric[, i], lwd = 2,
            ylim = c(0, 0.1),
            type = "l", 
            ylab = "Friction (unit)", xlab = "Domain (km)", 
            cex.axis = 3, cex.lab = 4)
        # lines(domain[plot_domain]/1000, true_fric[i, ], col = "red")
        lines(domain[plot_domain]/1000, pred_fric[, i], lwd = 2, col = "red")
        lines(domain[plot_domain]/1000, fric_lq[[i]], lwd = 2, lty = 1, col = "royalblue")
        lines(domain[plot_domain]/1000, fric_uq[[i]], lwd = 2, lty = 1, col = "royalblue")
        abline(v = test_gl[samples[i], years], lty = 2, lwd = 3)
    }

    dev.off()

    # ## Friction plot (diagonal)
    # png(paste0(plot_dir, "/diag_pred_vs_true_fric", checkpt, ".png"), width = 2000, height = 1200)

    # par(mfrow = c(length(samples) / 2, 2))
    # for (i in 1:length(samples)) {
    #     par(mar = c(6, 8, 4, 2))
    #     gl <- test_data$grounding_line[i] / 800 * 2001
    #     plot_domain <- 1:length(domain) # ceiling(gl)

    #     plot(test_fric[, i], pred_fric[, i], 
    #         ylab = "Predicted friction (unit)", xlab = "True friction (unit)", 
    #         cex.axis = 3, cex.lab = 4)
    #     # lines(domain[plot_domain]/1000, fric_lq[[i]], lty = 1, col = "royalblue")
    #     # lines(domain[plot_domain]/1000, fric_uq[[i]], lty = 1, col = "royalblue")
    #     abline(a = 0, b = 1, col = "red")
    # }

    # dev.off()

    ## Friction plot on real data
    png(paste0(plot_dir, "/pred_fric_real.png"), width = 1000, height = 500)
    plot(domain/1e3, pred_fric_real, type = "l", ylim = c(0, 0.1),
        ylab = "Friction (unit)", xlab = "Domain (km)", 
        cex.axis = 2, cex.lab = 2)
    lines(domain/1e3, fric_lq_real, lty = 1, col = "royalblue")
    lines(domain/1e3, fric_uq_real, lty = 1, col = "royalblue")
    dev.off()

    ## Bed plot
    png(paste0(plot_dir, "/pred_vs_true_bed.png"), width = 2000, height = 1200)
    par(mfrow = c(length(samples) / 2, 2))

    for (i in 1:length(samples)) {
        par(mar = c(6, 8, 4, 2))
        # gl <- test_data$grounding_line[i] / 800 * 2001
        plot_domain <- 1:length(domain) # ceiling(gl)
        plot(domain[plot_domain] / 1000, true_bed[i, plot_domain],
            type = "l", lwd = 2, col = "black",
            ylab = "Bed (m)", xlab = "Domain (km)", cex.axis = 3, cex.lab = 4
        )
        # lines(domain[plot_domain] / 1000, test_bed[plot_domain, i], lwd = 2, col = "blue")
        lines(domain[plot_domain] / 1000, pred_bed[plot_domain, i], lwd = 3, col = "red")
        lines(domain[plot_domain] / 1000, bed_lq[[i]], lty = 1, lwd = 3, col = "royalblue")
        lines(domain[plot_domain] / 1000, bed_uq[[i]], lty = 1, lwd = 3, col = "royalblue")
        points(bed_obs_chosen$loc/1000, bed_obs_chosen$bed_elev, pch = 16, cex = 1.5, col = "turquoise")
        abline(v = test_gl[samples[i], years], lty = 2, lwd = 3)
    }
    dev.off()

    # ## Bed plot (diagonal)
    # png(paste0(plot_dir, "/diag_pred_vs_true_bed", checkpt, ".png"), width = 2000, height = 1200)

    # par(mfrow = c(length(samples) / 2, 2))
    # for (i in 1:length(samples)) {
    #     par(mar = c(6, 8, 4, 2))
    #     gl <- test_data$grounding_line[i] / 800 * 2001
    #     plot_domain <- 1:length(domain) # ceiling(gl)

    #     plot(test_bed[, i], pred_bed[, i], 
    #         ylab = "Predicted bed elevation (m)", xlab = "True bed elevation (m)", 
    #         cex.axis = 3, cex.lab = 4)
    #     # lines(domain[plot_domain]/1000, fric_lq[[i]], lty = 1, col = "royalblue")
    #     # lines(domain[plot_domain]/1000, fric_uq[[i]], lty = 1, col = "royalblue")
    #     abline(a = 0, b = 1, col = "red")
    # }

    # dev.off()

    ## Bed plot on real data
    png(paste0(plot_dir, "/pred_bed_real.png"), width = 1000, height = 500)
    plot_domain <- 1:length(domain) #1500
    plot(domain[plot_domain]/1e3, pred_bed_real[plot_domain], type = "l", 
        ylab = "Bed elevation (m)", xlab = "Domain (km)", 
        cex.axis = 2, cex.lab = 2)
    lines(domain[plot_domain]/1e3, bed_lq_real[plot_domain], lty = 1, col = "salmon")
    lines(domain[plot_domain]/1e3, bed_uq_real[plot_domain], lty = 1, col = "salmon")
    points(bed_obs_chosen$loc/1e3, bed_obs_chosen$bed_elev, pch = 16, cex = 1.5, col = "turquoise")
    points(bed_obs_validate$loc/1e3, bed_obs_validate$bed_elev, pch = 1, cex = 1.5, col = "red")
    dev.off()

    ## Grounding line plot
    years <- dim(test_data$input)[3]
    png(paste0(plot_dir, "/pred_vs_true_gl.png"), width = 2000, height = 2000)
    par(mfrow = c(length(samples) / 2, 2))

    for (i in 1:length(samples)) {
        par(mar = c(6, 8, 4, 2))
        # gl <- test_data$grounding_line[i] / 800 * 2001
        plot_domain <- 1:length(domain) # ceiling(gl)
        plot(true_gl[i, ], 1:years,
            type = "l", lwd = 2, 
            xlab = "Grounding line (km)", ylab = "Time (year)", cex.axis = 3, cex.lab = 4
        )
        # lines(domain[plot_domain] / 1000, test_fric[i, ], lwd = 2, col = "blue")
        lines(pred_gl[i, ], 1:years, lwd = 2, col = "red")
        lines(pred_gl_lq[i, ], 1:years, lty = 1, lwd = 3, col = "royalblue")
        lines(pred_gl_uq[i, ], 1:years, lty = 1, lwd = 3, col = "royalblue")
        
        # abline(v = test_data$test_gl[i], lwd = 3, lty = 2)
    }

    dev.off()

    # png(paste0(plot_dir, "/diag_pred_vs_true_gl.png"), width = 2000, height = 1200)

    # par(mfrow = c(length(samples) / 2, 2))
    # for (i in 1:length(samples)) {
    #     par(mar = c(6, 8, 4, 2))
    #     # gl <- test_data$grounding_line[i] / 800 * 2001
    #     # plot_domain <- 1:length(domain) # ceiling(gl)

    #     plot(test_gl[, i], pred_gl[, i],  
    #         ylab = "Predicted bed elevation (m)", xlab = "True bed elevation (m)", 
    #         cex.axis = 3, cex.lab = 4)
    #     # lines(domain[plot_domain]/1000, fric_lq[[i]], lty = 1, col = "royalblue")
    #     # lines(domain[plot_domain]/1000, fric_uq[[i]], lty = 1, col = "royalblue")
    #     abline(a = 0, b = 1, col = "red")
    # }

    # dev.off()
}

############################################# RMSE ##################################################
rmse <- function(truth, pred) {
    sqrt(mean((truth - pred)^2))
}

total_rmse <- 0
# if (output_var == "friction") {## Post-process output from CNN
#     truth <- test_fric
#     pred <- pred_fric
# } else if (output_var == "grounding_line") {
#    truth <- test_gl
#    pred <- pred_gl
# } else if (output_var == "bed") {
#     truth <- test_bed
#     pred <- pred_bed
# } else if (output_var == "all") {
    truth <- cbind(test_fric, test_bed)
    pred <- cbind(pred_fric, pred_bed)
# } else {
#     print("Output variable should be either 'friction', 'bed' or 'grounding_line'.")
# }

for (s in 1:ncol(truth)) {
    total_rmse <- total_rmse + rmse(truth[, s], pred[, s])
}
cat("RMSE: ", total_rmse / ncol(truth), "\n") # normalise by number of simulations in the test set

# }
