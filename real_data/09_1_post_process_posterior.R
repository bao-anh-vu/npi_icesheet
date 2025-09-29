## Post-process output from CNN
setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())
source("./source/create_model.R")
source("./source/posterior_loss.R")
source("./source/sample_from_posterior.R")

library(parallel)
library(mvtnorm) 
library(Matrix)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(qs)
library(dplyr)

## Flags
save_pred <- T
save_plots <- T
log_transform <- T
test_on_train <- F
use_missing_pattern <- T

## Read data
data_date <- "20241111" #"20241103"
sets <- 51:100 #6:20
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Directories
# if (use_missing_pattern) {
#     data_dir <- paste0("./data/training_data/", setsf, "/")
#     output_dir <- paste0("./output/cnn/", setsf, "/")
#     plot_dir <- paste0("./plots/cnn/", setsf, "/")
# } else {
#     data_dir <- paste0("./data/training_data/", setsf, "/non")
#     output_dir <- paste0("./output/cnn/", setsf, "/non")
#     plot_dir <- paste0("./plots/cnn/", setsf, "/non")
# }
# if (correct_model_discrepancy) {
#     data_dir <- paste0("./data/training_data/", setsf, "/discr")
#     output_dir <- paste0("./output/cnn/", setsf, "/discr")
#     plot_dir <- paste0("./plots/cnn/", setsf, "/discr")
# } else {
    data_dir <- paste0("./data/training_data/", setsf, "/")
    output_dir <- paste0("./output/cnn/", setsf, "/")
    pred_output_dir <- paste0(output_dir, "/pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/")
# }

if (test_on_train) {
    test_data <- qread(file = paste0(data_dir, "/train_data_", data_date, ".qs"))
} else {
    test_data <- qread(file = paste0(data_dir, "/test_data_", data_date, ".qs"))
}   

test_input <- test_data$input
test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)

if (test_on_train) { ## Restrict input and output to the first 100 samples
    test_input <- test_input[1:100,,,]
    test_output <- test_output[1:100, ]
}

setf <- formatC(sets[1], width = 2, flag = "0")
friction_basis <- qread(file = paste0("./data/training_data/friction_basis_", setf, "_", data_date, ".qs"))
dim(friction_basis$basis_mat)      

# png(paste0(output_dir, "/true_fric_1.png"), width = 2000, height = 1200)
# plot(friction_basis$fitted_values[1, ], type = "l")
# dev.off()

ssa_steady <- qread(file = paste0("data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain

## Bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
bed_obs_train <- bed_obs_df %>% filter(chosen == 1) 
bed_obs_val <- bed_obs_df %>% filter(chosen == 0) 

bed_basis <- qread(file = paste0("./data/training_data/bed_basis_", setf, "_", data_date, ".qs"))
bed_mean <- bed_basis$mean 

## Load the model

input_dim <- dim(test_data$input)[2:4]
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

# ### Plot the loss
# history <- qread(file = paste0(output_dir, "/history_", data_date, ".qs"))
# # history %>% plot() #+
# # coord_cartesian(xlim = c(1, epochs))

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
history <- qread(file = paste0(output_dir, "history_", data_date, ".qs"))
checkpoint_path <- paste0(output_dir, "checkpoints/cp-{epoch:04d}.ckpt")
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
# if (test_on_train) {
#     train_subset <- train_data$input[1:100,,,]
#     pred_output <- model %>% predict(train_subset)
# } else {
    pred_time <- system.time({
        pred_output <- model %>% predict(test_input)
    })
    cat("Prediction time: ", pred_time[3], "\n")
# }

## Plot the predicted friction coefs generated by these basis functions
friction_basis <- qread(file = paste0("data/training_data/friction_basis_", setf, "_", data_date, ".qs"))
fric_basis_mat <- friction_basis$basis_mat
n_fric_basis <- ncol(fric_basis_mat)

bed_basis <- qread(file = paste0("data/training_data/bed_basis_", setf, "_", data_date, ".qs"))
bed_basis_mat <- bed_basis$basis_mat
n_bed_basis <- ncol(bed_basis_mat)

pred_mean <- pred_output[, 1:n_mean_elements]

## Un-standardise output
pred_fric_coefs <- pred_mean[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
pred_bed_coefs <- pred_mean[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
pred_gl <- pred_mean[, (n_fric_basis+n_bed_basis+1):ncol(pred_mean)] * test_data$sd_gl + test_data$mean_gl

# if (test_on_train) {
#     train_fric_coefs <- train_data$fric_coefs[1:100, ] * train_data$sd_fric_coefs + train_data$mean_fric_coefs
#     train_bed_coefs <- train_data$bed_coefs[1:100, ] * train_data$sd_bed_coefs + train_data$mean_bed_coefs
# }
test_fric_coefs <- test_output[, 1:n_fric_basis] * test_data$sd_fric_coefs + test_data$mean_fric_coefs
test_bed_coefs <- test_output[, (n_fric_basis+1):(n_fric_basis+n_bed_basis)] * test_data$sd_bed_coefs + test_data$mean_bed_coefs
test_gl <- test_output[, (n_fric_basis+n_bed_basis+1):ncol(test_output)] * test_data$sd_gl + test_data$mean_gl

if (test_on_train) {
    plot_tag <- "_train"
} else {
    plot_tag <- "_test"
}

## Plot predicted vs true basis coefficients for several samples
samples <- 1:5
png(paste0(plot_dir, "pred_coef", plot_tag, ".png"), width = 1000, height = 200*length(samples))
par(mfrow = c(length(samples), 2))

for (s in samples) {
    plot(test_fric_coefs[s,], 
        type = "l", main = "Friction basis coefficients")
    lines(pred_fric_coefs[s,], col = "red")

    plot(test_bed_coefs[s,], 
        type = "l", main = "Bed basis coefficients")
    lines(pred_bed_coefs[s,], col = "red")
}
dev.off()

## Compute predicted vs original friction 
if (log_transform) {
    pred_fric <- exp(fric_basis_mat %*% t(pred_fric_coefs))
    test_fric <- exp(fric_basis_mat %*% t(test_fric_coefs))
} else {
    pred_fric <- fric_basis_mat %*% t(pred_fric_coefs)
    test_fric <- fric_basis_mat %*% t(test_fric_coefs)
}

## Compute predicted vs original bed fluctuations (demeaned bed)
pred_bed_demean <- bed_basis_mat %*% t(pred_bed_coefs)
test_bed_demean <- bed_basis_mat %*% t(test_bed_coefs)

### Reconstruct bed elevation by adding bed trend to the predicted oscillations
# bed_obs <- readRDS(file = paste0("./data/training_data/bed_obs_", data_date, ".rds"))
# bed_obs <- readRDS(file = "./data/bed_elev_nearest.rds")
# bed_obs <- bed_obs[1:years01]

# bed_df <- data.frame(obs_locations = domain, bed_elev = bed_obs)
# bed_df <- na.omit(bed_df)
# bed.fit <- loess(bed_elev ~ loc, data = bed_obs, span = 0.25, 
#                     control = loess.control(surface = "direct")) 
# bed_mean <- predict(bed.fit, newdata = data.frame(loc = domain))
bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = ncol(pred_bed_demean))

pred_bed <- pred_bed_demean + bed_mean_mat
test_bed <- test_bed_demean + bed_mean_mat

###################################
##   Uncertainty quantification  ##
###################################
print("Computing credible intervals...")

## Covariance matrix
pred_chol <- pred_output[, (n_mean_elements+1):ncol(pred_output)]

## Construct Cholesky factor of the precision
Lb_elems <- n_bed_basis + (n_bed_basis - 1)
Lc_elems <- n_fric_basis + (n_fric_basis - 1)
Lg_elems <- n_gl + (n_gl - 1)

Lmats <- list()
s0 <- system.time({
    for (s in 1:nrow(pred_chol)) {
        Lb <- construct_L_matrix(pred_chol[s, 1:Lb_elems], n_bed_basis)
        Lc <- construct_L_matrix(pred_chol[s, (Lb_elems+1):(Lb_elems+Lc_elems)], n_fric_basis)
        Lg <- construct_L_matrix(pred_chol[s, (Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)
        Lmats[[s]] <- bdiag(Lb, Lc, Lg)
    }
})


# L1 <- Lmats[[1]]
# L1_inv <- solve(L1)
# Q <- t(L1_inv) %*% L1_inv
# image(Q)

# ## Need to sample from the posterior distribution of the coefs
# ## then transform them to actual friction, bed, gl

S <- 1000 ## number of posterior samples

## Sample basis function coefficients from posterior
s1 <- system.time({
pred_samples_ls <- lapply(1:nrow(test_output), function(i) 
                            sample_from_posterior(n = S, mean = pred_mean[i, ], prec_chol = Lmats[[i]])) #,
                                # mc.cores = 10L)
})  

## Transform basis coefficients into actual friction coefficients
sd_fric_coefs <- test_data$sd_fric_coefs
mean_fric_coefs <- test_data$mean_fric_coefs

s2 <- system.time({
    fric_samples_ls <- lapply(pred_samples_ls, function(x, log_transform, mean_fric_coefs, sd_fric_coefs) { 
    fric_coefs <- x[1:n_fric_basis, ]
    fric_coefs_ustd <- fric_coefs * sd_fric_coefs + mean_fric_coefs
    
    # if (log_transform) {
        fric_samples <- exp(fric_basis_mat %*% fric_coefs_ustd)
    # } else {
        # fric_samples <- fric_basis_mat %*% fric_coefs_ustd
    # }
    return(fric_samples)
    }, 
    log_transform = log_transform,
    mean_fric_coefs = mean_fric_coefs,
    sd_fric_coefs = sd_fric_coefs#,
    # mc.cores = 4L
)

})


## Transform basis coefficients into bed elevations
bed_mean_mat <- matrix(rep(bed_mean), nrow = length(bed_mean), ncol = S)
s3 <- system.time({
bed_samples_ls <- lapply(pred_samples_ls, function(x) {
    bed_coefs <- x[(n_fric_basis+1):(n_fric_basis+n_bed_basis), ]
    bed_coefs_ustd <- bed_coefs * test_data$sd_bed_coefs + test_data$mean_bed_coefs
    bed_samples <- bed_basis_mat %*% bed_coefs_ustd + bed_mean_mat
    }
)

})

## Grounding line is a direct output from the CNN
gl_samples_ls <- lapply(pred_samples_ls, function(x) x[(n_fric_basis+n_bed_basis+1):n_mean_elements, ])

## Compute quantiles
s5 <- system.time({
    fric_q <- mclapply(fric_samples_ls, function(x) apply(x, 1, quantile, probs = c(0.025, 0.975)), mc.cores = 10L)
    bed_q <- mclapply(bed_samples_ls, function(x) apply(x, 1, quantile, probs = c(0.025, 0.975)), mc.cores = 10L)
    gl_q <- mclapply(gl_samples_ls, function(x) apply(x, 1, quantile, probs = c(0.025, 0.975)), mc.cores = 10L)
})

fric_lq <- lapply(fric_q, function(x) x[1, ])
fric_uq <- lapply(fric_q, function(x) x[2, ])
bed_lq <- lapply(bed_q, function(x) x[1, ])
bed_uq <- lapply(bed_q, function(x) x[2, ])
gl_lq <- lapply(gl_q, function(x) x[1, ])
gl_uq <- lapply(gl_q, function(x) x[2, ])


#################################
##      Save predictions       ##
#################################
print("Saving predictions...")
if (save_pred) {
    qsave(pred_fric, file = paste0(pred_output_dir, "/pred_fric_", data_date, ".qs"))
    qsave(pred_bed, file = paste0(pred_output_dir, "/pred_bed_", data_date, ".qs"))
    qsave(pred_gl, file = paste0(pred_output_dir, "/pred_gl_", data_date, ".qs"))
    qsave(Lmats, file = paste0(pred_output_dir, "/Lmats_", data_date, ".qs"))
    qsave(fric_samples_ls, file = paste0(pred_output_dir, "/fric_post_samples_", data_date, ".qs"))
    qsave(bed_samples_ls, file = paste0(pred_output_dir, "/bed_post_samples_", data_date, ".qs"))
    qsave(gl_samples_ls, file = paste0(pred_output_dir, "/gl_post_samples_", data_date, ".qs"))
}

######################################
##      Comparison with truth       ##
######################################

## True parameter values for comparison
true_fric <- t(test_data$true_fric)
true_bed <- test_data$true_bed
true_gl <- test_data$grounding_line * test_data$sd_gl + test_data$mean_gl

## Scaling units for friction coefficients
# secpera <- 31556926
# fric_scale <- 1e6 * secpera^(1/3)
# true_fric <- true_fric/fric_scale

## Choose random samples from test set for comparison
samples <- 1:6 #sample(1:nrow(test_output), 6)

if (save_plots) {
    print("Saving plots...")
    ## Friction plots
    png(paste0(plot_dir, "/pred_vs_true_fric", plot_tag, ".png"), width = 2000, height = 1500, res = 100)

    plot_domain <- 1:1500 #length(domain) # ceiling(gl)

    par(mfrow = c(length(samples) / 2, 2))
    # for (i in 1:length(samples)) {
    for (s in samples) {
        # s <- samples[i]
        par(mar = c(6, 8, 4, 2))
        gl <- test_data$grounding_line[s] / 800 * 2001

        # plot(domain[plot_domain]/1000, test_fric[plot_domain, s], type = "l", ylim = c(0, 0.1),
        #     ylab = "Friction (unit)", xlab = "Domain (km)", lwd = 3,
        #     cex.axis = 3, cex.lab = 4,
        #     main = paste0("Sample ", s))
        plot(domain[plot_domain]/1000, test_fric[plot_domain, s], type = "l", #ylim = c(0, 0.06),
            ylab = "Friction (unit)", xlab = "Domain (km)", lwd = 3,
            cex.axis = 3, cex.lab = 3, cex.main = 3,
            main = paste0("Sample ", s))
        matlines(domain[plot_domain]/1000, fric_samples_ls[[s]][plot_domain, 1:3], 
                    # ylim = c(0, 0.06),
                    lty = 1, lwd = 2, col = adjustcolor("mediumpurple", alpha = 0.5))
        # lines(domain[plot_domain]/1000, true_fric[i, ], col = "red")
        lines(domain[plot_domain]/1000, pred_fric[plot_domain, s], lwd = 3, col = "red")
        lines(domain[plot_domain]/1000, fric_lq[[s]][plot_domain], lty = 1, lwd = 2, col = "salmon")
        lines(domain[plot_domain]/1000, fric_uq[[s]][plot_domain], lty = 1, lwd = 2, col = "salmon")
        
        abline(v = true_gl[s, ncol(true_gl)], lty = 2, lwd = 3)
    }

    dev.off()


    ## Bed plot
    png(paste0(plot_dir, "/pred_vs_true_bed", plot_tag, ".png"), width = 2000, height = 1500, res = 100)
    par(mfrow = c(length(samples) / 2, 2))

    # for (i in 1:length(samples)) {
    for (s in samples) {
        # s <- samples[i]
        par(mar = c(6, 8, 4, 2))
        # gl <- test_data$grounding_line[i] / 800 * 2001
        plot(domain[plot_domain] / 1000, test_bed[plot_domain, s],
            type = "l", lwd = 3, col = "black",
            ylab = "Bed (m)", xlab = "Domain (km)", 
            cex.axis = 3, cex.lab = 3, cex.main = 3,
            main = paste0("Sample ", s)
        )
        # lines(domain[plot_domain] / 1000, test_bed[plot_domain, i], lwd = 2, col = "blue")
        matlines(domain[plot_domain]/1000, bed_samples_ls[[s]][plot_domain, 1:3], 
                lty = 1, lwd = 2, col = adjustcolor("mediumpurple", alpha = 0.5))
        lines(domain[plot_domain] / 1000, pred_bed[plot_domain, s], lwd = 3, col = "red")
        lines(domain[plot_domain] / 1000, bed_lq[[s]][plot_domain], lty = 1, lwd = 2, col = "salmon")
        lines(domain[plot_domain] / 1000, bed_uq[[s]][plot_domain], lty = 1, lwd = 2, col = "salmon")
        points(bed_obs_train$loc / 1000, bed_obs_train$bed_elev, col = "black", pch = 16, cex = 2)
        points(bed_obs_val$loc / 1000, bed_obs_val$bed_elev, col = "cyan", pch = 16, cex = 2)
        abline(v = true_gl[s, ncol(true_gl)], lty = 2, lwd = 3)
    }
    dev.off()

    ## Grounding line plot
    years <- dim(true_gl)[2]
    png(paste0(plot_dir, "/pred_vs_true_gl", plot_tag, ".png"), width = 2000, height = 1500, res = 100)
    par(mfrow = c(length(samples) / 2, 2))
    
    # for (i in 1:length(samples)) {
    for (s in samples) {
        # s <- samples[i]
        par(mar = c(6, 8, 4, 2))
        # gl <- test_data$grounding_line[i] / 800 * 2001
        plot_domain <- 1:length(domain) # ceiling(gl)
        plot(true_gl[s, ], 1:years,
            type = "l", lwd = 3, 
            xlab = "Grounding line (km)", ylab = "Time (year)", cex.axis = 3, cex.lab = 4
        )
        matlines(gl_samples_ls[[s]][, 1:3], 1:years, lty = 1, lwd = 2, col = adjustcolor("mediumpurple", alpha = 0.5))
        # lines(domain[plot_domain] / 1000, test_fric[i, ], lwd = 2, col = "blue")
        lines(pred_gl[s, ], 1:years, lwd = 3, col = "red")
        lines(gl_lq[[s]], 1:years, lty = 1, lwd = 3, col = "salmon")
        lines(gl_uq[[s]], 1:years, lty = 1, lwd = 3, col = "salmon")
        
        # abline(v = test_data$test_gl[i], lwd = 3, lty = 2)
    }

    dev.off()

}
