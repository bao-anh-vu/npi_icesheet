## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/")

# source("sim_data.R")
source("./source/create_model.R")
source("./source/posterior_loss.R")

library(mvtnorm)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(gridExtra)
library(qs)

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

## Read data
data_date <- "20250108"
simulated_data <- qread(paste0("output/lgss_data_", data_date, ".qs"))
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
n_mean_elements <- 3 # phi, sigma_eta, sigma_eps
n_chol_elements <- 5 # diagonal + lower-diag elements
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
pred_time <- system.time({
    pred_output <- model %>% predict(test_input)
})
cat("Prediction time: ", pred_time[3], "\n")

pred_mean <- pred_output[, 1:n_mean_elements]

# Compare predicted mean vs true parameters for a chosen sample in the test set
phi_pred <- tanh(pred_mean[, 1]) # back-transform
sigma_eta_pred <- sqrt(exp(pred_mean[, 2]))
sigma_eps_pred <- sqrt(exp(pred_mean[, 3]))

# test_sample <- test_output[sample, ] # arctanh(phi), log(sigma_eta), log(sigma_eps)
phi_true <- tanh(test_output[, 1])
sigma_eta_true <- sqrt(exp(test_output[, 2]))
sigma_eps_true <- sqrt(exp(test_output[, 3]))

phi_df <- data.frame(param = "phi", pred = phi_pred, true = phi_true)
sigma_eta_df <- data.frame(param = "sigma[eta]", pred = sigma_eta_pred, true = sigma_eta_true)
sigma_eps_df <- data.frame(param = "sigma[epsilon]", pred = sigma_eps_pred, true = sigma_eps_true)
mean_pred_df <- rbind(phi_df, sigma_eta_df, sigma_eps_df)

mean_pred_plot <- ggplot(mean_pred_df, aes(x = true, y = pred)) + 
  geom_point(color = "salmon") + 
  geom_abline(intercept = 0, slope = 1, color = "black") +
  facet_wrap(~param, scales = "free", labeller = label_parsed) +
  theme(aspect.ratio = 1) +
  # guides(fill="none") +
  xlab("True parameters") +
  ylab("Predicted parameters") +
  # coord_fixed() +
  # scale_x_continuous(expand = c(0, 0)) + 
  # scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(text = element_text(size = 30))

png(paste0(output_dir, "/lgss_pred_vs_test_plot_", data_date, ".png"), width = 1200, height = 400)
print(mean_pred_plot)  
dev.off()

## Variance of the predicted parameters
construct_L_mat <- function(v, d) {
  L_mat <- diag(exp(v[1:d]))
  diag(L_mat[-1, -d]) <- v[(d + 1):length(v)]
  return(L_mat)
}

pred_chol <- pred_output[, (n_mean_elements + 1):output_dim]
pred_chol_ls <- lapply(1:nrow(pred_chol), function(i) construct_L_mat(pred_chol[i, ], n_mean_elements))
L_inv_ls <- lapply(pred_chol_ls, solve)

## For each test sample, simulate from the predicted posterior distribution
s <- 1
S <- 1000

phi_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)
sigma_eta_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)
sigma_eps_quantiles <- matrix(NA, nrow = nrow(test_output), ncol = 2)

for (s in 1:nrow(test_output)) {
    post_samples <- rmvnorm(S, pred_mean[s, ], t(L_inv_ls[[s]]) %*% L_inv_ls[[s]])

    ### Back-transform
    phi_post_samples <- tanh(post_samples[, 1])
    sigma_eta_post_samples <- sqrt(exp(post_samples[, 2]))
    sigma_eps_post_samples <- sqrt(exp(post_samples[, 3]))

    phi_quantiles[s, ] <- quantile(phi_post_samples, c(0.025, 0.975))
    sigma_eta_quantiles[s, ] <- quantile(sigma_eta_post_samples, c(0.025, 0.975))
    sigma_eps_quantiles[s, ] <- quantile(sigma_eps_post_samples, c(0.025, 0.975))

}

## Then maybe plot all of the samples point-wise
phi_df$lower <- phi_quantiles[, 1]
phi_df$upper <- phi_quantiles[, 2]
sigma_eta_df$lower <- sigma_eta_quantiles[, 1]
sigma_eta_df$upper <- sigma_eta_quantiles[, 2]
sigma_eps_df$lower <- sigma_eps_quantiles[, 1]
sigma_eps_df$upper <- sigma_eps_quantiles[, 2]

phi_dw_plot <- phi_df %>% ggplot() + 
    geom_errorbar(aes(x = 1:nrow(phi_df), ymin = lower, ymax = upper), width = 0, color = "grey") +
    geom_point(aes(x = 1:nrow(phi_df), y = pred), color = "black") + 
    geom_point(aes(x = 1:nrow(phi_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(phi)) +
    theme_bw() +
    theme(text = element_text(size = 20)) 

sigma_eta_plot <- sigma_eta_df %>% ggplot() +
    geom_errorbar(aes(x = 1:nrow(sigma_eta_df), ymin = lower, ymax = upper), width = 0, color = "grey") +
    geom_point(aes(x = 1:nrow(sigma_eta_df), y = pred), color = "grey") +
    geom_point(aes(x = 1:nrow(sigma_eta_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(sigma[eta])) +
    theme_bw() +
    theme(text = element_text(size = 20))  

sigma_eps_plot <- sigma_eps_df %>% ggplot() +
    geom_errorbar(aes(x = 1:nrow(sigma_eps_df), ymin = lower, ymax = upper), width = 0, color = "grey") +
    geom_point(aes(x = 1:nrow(sigma_eps_df), y = pred), color = "grey") +
    geom_point(aes(x = 1:nrow(sigma_eps_df), y = true), color = "salmon") +
    xlab("Test sample") +
    ylab(bquote(sigma[epsilon])) +
    theme_bw() +
    theme(text = element_text(size = 20))

png(paste0(output_dir, "/lgss_quantile_plot_", data_date, ".png"), width = 900, height = 500)
grid.arrange(phi_dw_plot, sigma_eta_plot, sigma_eps_plot, nrow = 3)
dev.off()

## Save prediction output
pred_mean_ls <- list(phi = phi_pred, sigma_eta = sigma_eta_pred, sigma_eps = sigma_eps_pred)
qsave(pred_mean_ls, paste0(output_dir, "/lgss_pred_mean_", data_date, ".qs"))
qsave(pred_chol_ls, paste0(output_dir, "/lgss_pred_chol_", data_date, ".qs"))


# par(mfrow = c(3,1))
# hist(phi_post_samples, main = "phi")
# abline(v = phi_true[s], col = "black")
# abline(v = phi_pred[s], col = "red")

# hist(sigma_eta_post_samples, main = "sigma_eta")
# abline(v = sigma_eta_true[s], col = "black")
# abline(v = sigma_eta_pred[s], col = "red")

# hist(sigma_eps_post_samples, main = "sigma_eps")
# abline(v = sigma_eps_true[s], col = "black")
# abline(v = sigma_eps_pred[s], col = "red")
# dev.off()
# ## Predict on the ``real'' data
# y <- qread(paste0("output/y_", data_date, ".qs"))
# y_std <- (y - mean_input) / sd_input
# y_std <- array_reshape(y_std, c(1, length(y_std), 1))

# real_pred <- model %>% predict(y_std)
# real_pred_mean <- real_pred[, 1:n_mean_elements]

# phi_pred <- tanh(real_pred_mean[1])
# sigma_eta_pred <- exp(real_pred_mean[2])
# sigma_eps_pred <- exp(real_pred_mean[3])
