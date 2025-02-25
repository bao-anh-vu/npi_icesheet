## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/full/")

rm(list = ls())

library(mvtnorm)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(gridExtra)
library(qs)

# source("sim_data.R")
source("./source/create_model.R")
source("./source/posterior_loss.R")


#List physical devices
gpus <- tf$config$experimental$list_physical_devices('GPU')

if (length(gpus) > 0) {
  tryCatch({
    # Restrict TensorFlofw to only allocate 4GB of memory on the first GPU
    tf$config$experimental$set_virtual_device_configuration(
      gpus[[1]],
      list(tf$config$experimental$VirtualDeviceConfiguration(memory_limit=4096*10))
    )
    
    logical_gpus <- tf$config$experimental$list_logical_devices('GPU')
    
    print(paste0(length(gpus), " Physical GPUs,", length(logical_gpus), " Logical GPUs"))
  }, error = function(e) {
    # Virtual devices must be set before GPUs have been initialized
    print(e)
  })
}

## Flags
use_arctanh <- T # use either arctanh or logit
test_on_train <- F

## Read data
data_date <- "20250108"
simulated_data <- qread(paste0("data/lgss_data_", data_date, ".qs"))
train_input <- simulated_data$train_input
val_input <- simulated_data$val_input
test_input <- simulated_data$test_input
train_output <- simulated_data$train_output
val_output <- simulated_data$val_output
test_output <- simulated_data$test_output

## Create CNN model
print("Creating model...")
# Create a basic model instance
input_dim <- c(dim(train_input)[-1])
n_mean_elements <- 3L # phi, sigma_eta, sigma_eps
n_chol_elements <- n_mean_elements + (n_mean_elements-1) # diagonal + first sub-diag elements
output_dim <- n_mean_elements + n_chol_elements  

model <- create_model_posterior(input_dim = input_dim, 
                                output_dim = output_dim, 
                                d = n_mean_elements)

# Display the model's architecture
summary(model)

output_dir <- "./output/"

history <- qread(paste0(output_dir, "/history_", data_date, ".qs"))
# checkpoint_path <- paste0(output_dir, "/checkpoints/cp-{epoch:04d}.ckpt")

## Choose model with lowest validation loss
checkpt <- which.min(history$metrics$val_loss) 
cp_restart <- paste0(output_dir, "checkpoints/cp-", formatC(checkpt, width = 4, flag = "0"), ".ckpt")
# latest <- tf$train$latest_checkpoint(checkpoint_dir)
load_model_weights_tf(model, cp_restart)

## Predict on test set
if (test_on_train) {
    test_input <- train_input[1:100, , ]
    test_output <- train_output[1:100]
}
pred_time <- system.time({
    pred_output <- model %>% predict(test_input)
})
cat("Prediction time: ", pred_time[3], "\n")

pred_mean <- pred_output[, 1:n_mean_elements]
pred_chol <- pred_output[, (n_mean_elements + 1):output_dim]

## Variance of the predicted parameters
construct_L_mat <- function(v, d) {
  L_mat <- diag(exp(v[1:d]))
  diag(L_mat[-1, -d]) <- v[(d + 1):length(v)]
  return(L_mat)
}
pred_chol_ls <- lapply(1:nrow(pred_output), function(i) construct_L_mat(pred_chol[i, ], d = n_mean_elements))

L_inv_ls <- lapply(pred_chol_ls, solve)



# png(paste0("./plots/lgss_pred_vs_test_plot_", data_date, ".png"), width = 500, height = 500)
# plot(test_output, pred_mean, xlab = "True phi", ylab = "Predicted phi", col = "salmon")
# abline(a = 0, b = 1, col = "black")
# dev.off()


## For each test sample, simulate from the predicted posterior distribution
S <- 10000

phi_mean <- c()
phi_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)
phi_post_samples_ls <- list()

sigma_eta_mean <- c()
sigma_eta_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)
sigma_eta_post_samples_ls <- list()

sigma_eps_mean <- c()
sigma_eps_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)
sigma_eps_post_samples_ls <- list()

for (s in 1:nrow(test_output)) {
    post_samples <- rmvnorm(S, pred_mean[s, ], t(L_inv_ls[[s]]) %*% L_inv_ls[[s]])
   
    ### Back-transform
    if (use_arctanh) {
        phi_post_samples <- tanh(post_samples[, 1])
    } else { # use logit
        phi_post_samples <- exp(post_samples[, 1]) / (1 + exp(-post_samples[, 1]))
    }

    sigma_eta_samples <- sqrt(exp(post_samples[, 2]))
    sigma_eps_samples <- sqrt(exp(post_samples[, 3]))

    phi_mean[s] <- mean(phi_post_samples)
    phi_quantiles[s, ] <- quantile(phi_post_samples, c(0.025, 0.975))
    phi_post_samples_ls[[s]] <- phi_post_samples

    sigma_eta_mean[s] <- mean(sigma_eta_samples)
    sigma_eta_quantiles[s, ] <- quantile(sigma_eta_samples, c(0.025, 0.975))
    sigma_eta_post_samples_ls[[s]] <- sigma_eta_samples

    sigma_eps_mean[s] <- mean(sigma_eps_samples)
    sigma_eps_quantiles[s, ] <- quantile(sigma_eps_samples, c(0.025, 0.975))
    sigma_eps_post_samples_ls[[s]] <- sigma_eps_samples
}

if (use_arctanh) {
    phi_true <- tanh(test_output[, 1])
} else { # use logit
    phi_true <- exp(test_output[, 1]) / (1 + exp(-test_output[, 1]))
}
sigma_eta_true <- sqrt(exp(test_output[, 2]))
sigma_eps_true <- sqrt(exp(test_output[, 3]))

## Then maybe plot all of the samples point-wise
phi_df <- data.frame(pred = phi_mean, true = phi_true,
                     lower = phi_quantiles[, 1], upper = phi_quantiles[, 2])
sigma_eta_df <- data.frame(pred = sigma_eta_mean, true = sigma_eta_true,
                           lower = sigma_eta_quantiles[, 1], upper = sigma_eta_quantiles[, 2])
sigma_eps_df <- data.frame(pred = sigma_eps_mean, true = sigma_eps_true,
                           lower = sigma_eps_quantiles[, 1], upper = sigma_eps_quantiles[, 2])

## Plot the mean against true values
phi_mean_plot <- phi_df %>% ggplot() +
              geom_point(aes(x = true, y = pred), color = "salmon") +
              geom_abline(intercept = 0, slope = 1, color = "black") +
              xlab("True phi") +
              ylab("Predicted phi") +
              theme_bw() +
              theme(text = element_text(size = 20))

sigma_eta_mean_plot <- sigma_eta_df %>% ggplot() +
              geom_point(aes(x = true, y = pred), color = "salmon") +
              geom_abline(intercept = 0, slope = 1, color = "black") +
              xlab("True sigma_eta") +
              ylab("Predicted sigma_eta") +
              theme_bw() +
              theme(text = element_text(size = 20))


sigma_eps_mean_plot <- sigma_eps_df %>% ggplot() +
              geom_point(aes(x = true, y = pred), color = "salmon") +
              geom_abline(intercept = 0, slope = 1, color = "black") +
              xlab("True sigma_eps") +
              ylab("Predicted sigma_eps") +
              theme_bw() +
              theme(text = element_text(size = 20))

## Plot the mean against the true values
png(paste0("./plots/lgss_mean_plot_", data_date, ".png"), width = 1000, height = 500)
# print(mean_plot)
grid.arrange(phi_mean_plot, sigma_eta_mean_plot, sigma_eps_mean_plot, ncol = 2)
dev.off()

phi_intv_plot <- phi_df %>% ggplot() + 
    geom_errorbar(aes(x = 1:nrow(phi_df), ymin = lower, ymax = upper), width = 0, color = "black") +
    geom_point(aes(x = 1:nrow(phi_df), y = pred), color = "black") + 
    geom_point(aes(x = 1:nrow(phi_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(phi)) +
    theme_bw() +
    theme(text = element_text(size = 20)) 

sigma_eta_intv_plot <- sigma_eta_df %>% ggplot() + 
    geom_errorbar(aes(x = 1:nrow(sigma_eta_df), ymin = lower, ymax = upper), width = 0, color = "black") +
    geom_point(aes(x = 1:nrow(sigma_eta_df), y = pred), color = "black") + 
    geom_point(aes(x = 1:nrow(sigma_eta_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(sigma[epsilon])) +
    theme_bw() +
    theme(text = element_text(size = 20))

sigma_eps_intv_plot <- sigma_eps_df %>% ggplot() + 
    geom_errorbar(aes(x = 1:nrow(sigma_eps_df), ymin = lower, ymax = upper), width = 0, color = "black") +
    geom_point(aes(x = 1:nrow(sigma_eps_df), y = pred), color = "black") + 
    geom_point(aes(x = 1:nrow(sigma_eps_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(sigma[epsilon])) +
    theme_bw() +
    theme(text = element_text(size = 20))

png(paste0("./plots/lgss_intv_plot_", data_date, ".png"), width = 900, height = 500)
grid.arrange(phi_intv_plot, sigma_eta_intv_plot, sigma_eps_intv_plot, nrow = 2)
dev.off()


## Compare to output from HMC
test_sample <- 3:7
hmc.phi_samples <- list()
hmc.sigma_eta_samples <- list()
hmc.sigma_eps_samples <- list()

for (s in test_sample) {
  hmc_output <- qread(paste0("output/lgss_hmc_results_", 
                            formatC(s, width = 2, format = "d", flag = "0"), 
                            "_", data_date, ".qs"))
  hmc.phi_samples[[s]] <- hmc_output$draws[, , 1]
  hmc.sigma_eta_samples[[s]] <- hmc_output$draws[, , 2]
  hmc.sigma_eps_samples[[s]] <- hmc_output$draws[, , 3]
}

png(paste0("./plots/lgss_hmc_vs_cnn_", data_date, ".png"), width = 600, height = 1200)
par(mfrow = c(length(test_sample), 3))

for (s in test_sample) {
  plot(density(phi_post_samples_ls[[s]], na.rm = T), col = "red", main = "phi", xlab = "phi")
  lines(density(hmc.phi_samples[[s]], na.rm = T))
  abline(v = phi_true[s], col = "black", lty = 2)

  plot(density(sigma_eta_post_samples_ls[[s]], na.rm = T),  col = "red", main = "sigma_eta", xlab = "sigma_eta")
  lines(density(hmc.sigma_eta_samples[[s]], na.rm = T))
  abline(v = sigma_eta_true[s], col = "black", lty = 2)

  plot(density(sigma_eps_post_samples_ls[[s]], na.rm = T),  col = "red", main = "sigma_eps", xlab = "sigma_eps")
  lines(density(hmc.sigma_eps_samples[[s]], na.rm = T))
  abline(v = sigma_eps_true[s], col = "black", lty = 2)
}

dev.off()






