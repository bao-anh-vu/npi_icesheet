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
library(gridExtra)

## Flags
# save_pred <- T
save_plots <- T
log_transform <- T
test_on_train <- F
# use_missing_pattern <- T
correct_model_discrepancy <- T
correct_velocity_discrepancy <- T
avg_over_time <- F
leave_one_out <- T

## Read data
data_date <- "20241111" #"20241103"
sets <- 1:50 #51:100 #51:100 #6:20
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./data/training_data/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")

## Save predictions
if (correct_model_discrepancy) {
    if (avg_over_time) {
        pred_output_dir <- paste0(output_dir, "pred/discr_avg/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr_avg/")
    } else {
        pred_output_dir <- paste0(output_dir, "pred/discr/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")
    }
    
} else {
    pred_output_dir <- paste0(output_dir, "pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

if (!dir.exists(pred_output_dir)) {
    dir.create(pred_output_dir)
}

## Domain data
ssa_steady <- qread(file = paste0("data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
gl_pos <- qread(file = paste0("data/grounding_line/gl_pos.qs"))
gl_obs <- domain[gl_pos$ind]

## BedMachine data for comparison
bedmachine <- qread("./data/bedmachine/flowline_bedmachine.qs")

## Load test data
train_data <- qread(file = paste0(data_dir, "train_data_", data_date, ".qs"))
test_data <- qread(file = paste0(data_dir, "test_data_", data_date, ".qs"))
test_input <- test_data$input

n_years <- dim(test_input)[3]
if (leave_one_out) {
    n_years <- n_years - 1
    test_input <- test_input[,,1:n_years,]
}

test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)

## Load real surface observations (just use the first 10 years, leave the last for OOS forecast)
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

surf_elev_data <- surf_elev_data[, 1:n_years]
velocity_data <- velocity_data[, 1:n_years]

## Adjust real data if correcting for model discrepancy

## Model discrepancy
vel_discr <- qread(file = paste0("./data/discrepancy/vel_discr_", data_date, ".qs"))
se_discr <- qread(file = paste0("./data/discrepancy/se_discr_", data_date, ".qs"))

avg_vel_discr <- rowMeans(vel_discr, na.rm = T)
avg_se_discr <- rowMeans(se_discr, na.rm = T)
if (avg_over_time) { # just use average discrepancy for all years
    vel_discr_mat <- matrix(rep(avg_vel_discr, n_years), nrow = length(avg_vel_discr), ncol = n_years)
    se_discr_mat <- matrix(rep(avg_se_discr, n_years), nrow = length(avg_se_discr), ncol = n_years)
} else { 
    vel_discr_mat <- vel_discr
    se_discr_mat <- se_discr
}

if (correct_model_discrepancy) {
    surf_elev_data <- surf_elev_data - se_discr_mat
    if (correct_velocity_discrepancy) {
        velocity_data <- velocity_data - vel_discr_mat
    }  
}

## Replace NA values in real_data with 0
surf_elev_data[is.na(surf_elev_data)] <- 0
velocity_data[is.na(velocity_data)] <- 0

## Apply masked pattern to real data
print("Reading missing patterns...")
surf_elev_missing_pattern <- qread("./data/surface_elev/missing_pattern.qs")
vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")
# missing_patterns <- abind(list(surf_elev_missing_pattern, vel_missing_pattern), along = 3)

# if (leave_one_out) {
    surf_elev_missing_pattern <- surf_elev_missing_pattern[, 1:n_years]
    vel_missing_pattern <- vel_missing_pattern[, 1:n_years]
# }
surf_elev_data2 <- surf_elev_data * surf_elev_missing_pattern
velocity_data2 <- velocity_data * vel_missing_pattern

## Standardise real data using the mean and sd of the training/test data
input_mean <- test_data$input_mean
input_sd <- test_data$input_sd
mean_surf_elev <- input_mean[1] #mean(surf_elev_data, na.rm = T)
sd_surf_elev <- input_sd[1] #sd(surf_elev_data, na.rm = T)
surf_elev_data_std <- (surf_elev_data2 - mean_surf_elev) / sd_surf_elev

mean_velocity <- input_mean[2] #mean(velocity_data, na.rm = T)
sd_velocity <- input_sd[2] #sd(velocity_data, na.rm = T)
velocity_data_std <- (velocity_data2 - mean_velocity) / sd_velocity

# surf_elev_data_std2 <- surf_elev_data_std * surf_elev_missing_pattern
# velocity_data_std2 <- velocity_data_std * vel_missing_pattern
real_data <- abind(surf_elev_data_std, velocity_data_std, along = 3)
real_data <- array(real_data, dim = c(1L, dim(real_data)))

if (leave_one_out) {
    real_data <- real_data[,,1:n_years,]
    real_data <- array(real_data, dim = c(1L, dim(real_data))) # to make it match the dimension of test data
    print(paste0("Using first ", n_years, " years of data for leave-one-out cross-validation"))
}

## Compare real data and test data
# pdf("./plots/cnn/input/real_vs_test_data.pdf", width = 10, height = 10)
# par(mfrow = c(2,2))

# image(real_data[1,,,1], main = "Real surface elevation data")
# # axis(1, at = seq(0, 1, length.out = 5), labels = round(seq(min(domain)/1e3, max(domain)/1e3, length.out = 5),1), cex.axis = 1.5)
# # axis(2, at = seq(0, 1, length.out = 5), labels = round(seq(100, 1, length.out = 5),1), cex.axis = 1.5)
# mtext("Distance (km)", side = 1, line = 3, cex = 1.5)
# mtext("Time (yr)", side = 2, line = 3, cex = 1.5)

# image(test_input[1,,,1], main = "Test surface elevation data") #, col = terrain.colors(100))
# # axis(1, at = seq(0, 1, length.out = 5), labels = round(seq(min(domain)/1e3, max(domain)/1e3, length.out = 5),1), cex.axis = 1.5)
# # axis(2, at = seq(0, 1, length.out = 5), labels = round(seq(100, 1, length.out = 5),1), cex.axis = 1.5)  
# mtext("Distance (km)", side = 1, line = 3, cex = 1.5)
# mtext("Time (yr)", side = 2, line = 3, cex = 1.5)

# image(real_data[1,,,2], main = "Real velocity data")
# image(test_input[1,,,2], main = "Test velocity data") #, col = terrain.colors(100))
# dev.off()

## Same plot but using ggplot
space <- domain / 1000
time <- 1:dim(real_data)[3]
grid <- expand.grid(space, time)
head(grid)
names(grid) <- c("space", "time")

sim <- sample(1:dim(test_input)[1], 1)
grid$real_se <- as.vector(real_data[1,,,1])
grid$test_se <- as.vector(test_input[sim,,,1])
grid$real_vel <- as.vector(real_data[1,,,2])
grid$test_vel <- as.vector(test_input[sim,,,2])

se_lims <- range(c(grid$real_se, grid$test_se), na.rm = T)
vel_lims <- range(c(grid$real_vel, grid$test_vel), na.rm = T)

real_se_plot <- ggplot(grid) +
    geom_tile(aes(space, time, fill = real_se)) +
    scale_y_reverse() +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = se_lims) +
    theme_bw() +
    theme(text = element_text(size = 30)) +
    # labs(fill="Thickness (m)")
    labs(fill = "Surface elev. (m)") +
    ggtitle("Real data")

real_vel_plot <- ggplot(grid) +
    geom_tile(aes(space, time, fill = real_vel)) +
    scale_y_reverse() +
    theme_bw() +
    theme(text = element_text(size = 30)) +
    scale_fill_distiller(palette = "Reds", direction = 1, limits = vel_lims) +
    labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")")) +
    ggtitle("Real data")

test_se_plot <- ggplot(grid) +
    geom_tile(aes(space, time, fill = test_se)) +
    scale_y_reverse() +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = se_lims) +
    theme_bw() +
    theme(text = element_text(size = 30)) +
    # labs(fill="Thickness (m)")
    labs(fill = "Surface elev. (m)") +
    ggtitle("Test data")

test_vel_plot <- ggplot(grid) +
    geom_tile(aes(space, time, fill = test_vel)) +
    scale_y_reverse() +
    theme_bw() +
    theme(text = element_text(size = 30)) +
    scale_fill_distiller(palette = "Reds", direction = 1, limits = vel_lims) +
    labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")")) +
    ggtitle("Test data")

real_vs_test_plots <- list(test_se_plot, real_se_plot, test_vel_plot, real_vel_plot)

png(paste0(plot_dir, "real_vs_test_hovmoller.png"), width = 1400, height = 800)
par(mfrow = c(2,2))
grid.arrange(grobs = real_vs_test_plots)
dev.off()

## Compare real data and test data for some simulations
png(paste0(plot_dir, "real_vs_test_obs.png"), width = 800, height = 800, res = 100)
nsims <- 4
test_sims <- sample(1:dim(test_input)[1], nsims)
par(mfrow = c(nsims, 2))

for (sim in test_sims) {
    matplot(test_input[sim,,,1], col = "grey", type = "l", main = paste0("Simulation ", sim))
    matlines(real_data[1,,,1], col = "salmon")
    legend("topright", legend = c("test data", "real data"), col = c("grey", "salmon"), lty = 1)

    matplot(test_input[sim,,,2], col = "grey", type = "l", main = paste0("Simulation ", sim))
    matlines(real_data[1,,,2], col = "salmon")
    legend("topright", legend = c("test data", "real data"), col = c("grey", "salmon"), lty = 1)
}
dev.off()

## Load the model
input_dim <- dim(test_input)[2:4]
n_fric_basis <- dim(test_data$fric_coefs)[2]
n_bed_basis <- dim(test_data$bed_coefs)[2]
n_gl <- dim(test_data$grounding_line)[2]
n_mean_elements <- n_fric_basis + n_bed_basis + n_gl
n_chol_elements <- n_mean_elements + (n_mean_elements - 3) # diagonal + lower-diag elements
output_dim <- n_mean_elements + n_chol_elements  # THIS NEEDS TO CHANGE TO THE TOTAL NUMBER OF BASIS FUNCTIONS + COVARIANCE PARAMETERS

model <- create_model_posterior(input_dim = input_dim, 
                                output_dim = output_dim,
                                n_bed_basis = n_bed_basis,
                                n_fric_basis = n_fric_basis,
                                n_gl = n_gl)

### Display the model's architecture
summary(model)

### Plot the loss
history <- qread(file = paste0(output_dir, "history_", data_date, ".qs"))

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
    train_subset <- train_data$input[1,,1:10,]
    train_subset <- array(train_subset, dim = c(1L, dim(train_subset)))
    pred_output <- model %>% predict(train_subset)
} else {

    pred_time <- system.time({
    # pred_output <- model %>% predict(test_input)
    pred_output <- model %>% predict(real_data)

    })
    cat("Prediction time: ", pred_time[3], "\n")
}

## Plot the predicted friction coefs generated by these basis functions

# fric_basis_mat <- test_data$fric_basis_mat
# bed_basis_mat <- test_data$bed_basis_mat
# n_fric_basis <- ncol(fric_basis_mat)
# n_bed_basis <- ncol(bed_basis_mat)
setf <- formatC(sets[1], width=2, flag="0")
friction_basis <- qread(file = paste0("data/training_data/friction_basis_", setf, "_", data_date, ".qs"))
fric_basis_mat <- friction_basis$basis_mat
n_fric_basis <- ncol(fric_basis_mat)


bed_basis <- qread(file = paste0("data/training_data/bed_basis_", setf, "_", data_date, ".qs"))
bed_basis_mat <- bed_basis$basis_mat
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

png(paste0(plot_dir, "pred_coefs_real.png"), width = 2000, height = 1200)
par(mfrow = c(2,1))
plot(pred_fric_coefs, type = "l")

# png(paste0(plot_dir, "pred_bed_coefs_real.png"), width = 2000, height = 1200)
plot(pred_bed_coefs, type = "l")
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

# bed_basis <- qread(file = paste0("./data/training_data/bed_basis_01_", data_date, ".qs"))

bed_mean <- bed_basis$mean
# bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = ncol(pred_bed_demean))
pred_bed <- pred_bed_demean + bed_mean
# test_bed <- test_bed_demean + bed_mean_mat

## Covariance matrix
pred_chol <- pred_output[, (n_mean_elements+1):ncol(pred_output)]

## Construct Cholesky factor of the precision
Lb_elems <- n_bed_basis + (n_bed_basis - 1)
Lc_elems <- n_fric_basis + (n_fric_basis - 1)
Lg_elems <- n_gl + (n_gl - 1)

Lmats <- list()
# for (s in 1:nrow(pred_chol)) {
    Lb <- construct_L_matrix(pred_chol[1:Lb_elems], n_bed_basis)
    Lc <- construct_L_matrix(pred_chol[(Lb_elems+1):(Lb_elems+Lc_elems)], n_fric_basis)
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
bed_mean_mat <- matrix(rep(bed_mean, S), nrow = length(bed_mean), ncol = S)
bed_samples <- bed_basis_mat %*% bed_coefs_ustd + bed_mean_mat

## Grounding line is a direct output from the CNN
gl_samples_ls <- pred_samples_ls[(n_fric_basis+n_bed_basis+1):n_mean_elements, ]
gl_ustd <- gl_samples_ls * test_data$sd_gl + test_data$mean_gl

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
# png(paste0(plot_dir, "pred_fric_real.png"), width = 1000, height = 500)
# plot(domain/1e3, pred_fric, type = "l", col = "red", lwd = 2,
#     ylim = c(0, 0.1),
#      xlab = "Flowline (km)", ylab = "Friction coefficient")
# lines(domain/1e3, fric_uq, col = "grey", lwd = 2)
# lines(domain/1e3, fric_lq, col = "grey", lwd = 2)
# abline(v = gl_obs/1e3, lty = 2)
# dev.off()

## Same plot but using ggplot
fric_df <- data.frame(domain = domain/1e3, 
                      lq = fric_lq, 
                      uq = fric_uq, 
                      pred = pred_fric)
gl_df <- data.frame(gl = gl_obs/1e3)
fric_p <- ggplot(data = fric_df, aes(x = domain)) +
    geom_ribbon(aes(ymin = lq, ymax = uq), fill = "salmon", alpha = 0.5) +
    geom_line(aes(y = pred), color = "red", lwd = 1) +
    geom_vline(data = gl_df, aes(xintercept = gl), linetype = "dashed") +
    xlim(0, 150) +
    # ylim(0, 0.5) +
    labs(x = "Flowline (km)", y = expression(paste("Friction coefficient (", M ~ Pa ~ m^{-1/3} ~ a^{1/3}, ")"))) +
    theme_bw() +
    theme(text = element_text(size = 20),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))

# png(filename = paste0(plot_dir, "pred_fric_real_ggplot.png"), width = 1000, height = 500, res = 100)
# print(fric_p)
# dev.off()

## Validate against the rest of the bed obs
bed_obs_df <- qread(file = "./data/bedmap/bed_obs_df_all.qs")
# bed_obs_train <- bed_obs_df %>% filter(chosen == 1)
# bed_obs_val <- bed_obs_df %>% filter(chosen == 0)

# png(paste0(plot_dir, "pred_bed_real.png"), width = 1000, height = 500)
# plot(domain/1e3, bed_lq, type = "l", col = "grey", lwd = 2, xlab = "Flowline (km)", ylab = "Elevation (m)")
# lines(domain/1e3, bed_uq, col = "grey", lwd = 2)
# lines(domain/1e3, pred_bed, col = "red", lwd = 2)
# points(bed_obs_train$loc/1e3, bed_obs_train$bed_elev, col = "black", pch = 20)
# points(bed_obs_val$loc/1e3, bed_obs_val$bed_elev, col = "cyan")
# abline(v = gl_obs/1e3, lty = 2)
# dev.off()

## Same plot but using ggplot
bed_df <- data.frame(domain = domain/1e3, 
                     lq = bed_lq, 
                     uq = bed_uq, 
                     pred = pred_bed,
                     bedmachine = bedmachine$bed_avg,
                     prior_mean = bed_mean)
gl_df <- data.frame(gl = gl_obs/1e3)

bed_p <- ggplot(data = bed_df, aes(x = domain)) +
    geom_ribbon(aes(ymin = lq, ymax = uq), fill = "salmon", alpha = 0.5) +
    geom_line(aes(y = pred, color = "Predicted bed"), lwd = 1) +
    geom_line(aes(y = prior_mean, color = "Prior mean"), lwd = 1, lty = 3) +
    geom_line(aes(y = bedmachine, color = "BedMachine"), lwd = 1, lty = 2) +
    geom_point(data = bed_obs_df, aes(x = loc/1e3, y = bed_elev, color = "Observations"), size = 2) +
    geom_vline(data = gl_df, aes(xintercept = gl), linetype = "dashed") +
    xlim(0, 150) +
    ylim(-1500, -500) +
    labs(x = "Flowline (km)", y = "Elevation (m)", color = "") +
    scale_color_manual(
        values = c(
            "Predicted bed" = "red",
            "Prior mean" = "black",
            "BedMachine" = "blue",
            "Observations" = "black"
        )
    ) +
    theme_bw() +
    theme(
        text = element_text(size = 20),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        legend.position = "none"   # or "bottom"
    )


# ggsave(filename = paste0(plot_dir, "pred_bed_real_ggplot.png"), plot = p, width = 10, height = 5)
png(filename = paste0(plot_dir, "pred_bed_fric_real_ggplot.png"), width = 2000, height = 600, res = 120)
grid.arrange(fric_p, bed_p, ncol = 2)
dev.off()


# if (save_pred) {
    qsave(pred_fric, file = paste0(pred_output_dir, "pred_fric_real_", data_date, ".qs"))
    qsave(pred_bed, file = paste0(pred_output_dir, "pred_bed_real_", data_date, ".qs"))
    qsave(fric_samples, file = paste0(pred_output_dir, "/fric_samples_real_", data_date, ".qs"))
    qsave(bed_samples, file = paste0(pred_output_dir, "/bed_samples_real_", data_date, ".qs"))
    qsave(fric_q, file = paste0(pred_output_dir, "/fric_quantiles_real_", data_date, ".qs"))
    qsave(bed_q, file = paste0(pred_output_dir, "/bed_quantiles_real_", data_date, ".qs"))
    qsave(gl_q, file = paste0(pred_output_dir, "/gl_quantiles_real_", data_date, ".qs"))
# }