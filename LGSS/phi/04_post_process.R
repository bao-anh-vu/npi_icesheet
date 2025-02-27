## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/phi/")

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
library(dplyr)

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
output_dim <- 2 # mean and sd

model <- create_model_posterior(input_dim = input_dim, 
                                output_dim = output_dim)

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

pred_mean <- pred_output[, 1]
pred_sd <- exp(pred_output[, 2])

png(paste0("./plots/lgss_pred_vs_test_plot_", data_date, ".png"), width = 500, height = 500)
plot(test_output, pred_mean, xlab = "True phi", ylab = "Predicted phi", col = "salmon")
abline(a = 0, b = 1, col = "black")
dev.off()


# Compare predicted mean vs true parameters for a chosen sample in the test set
# phi_pred <- tanh(pred_mean[, 1]) # back-transform
# sigma_eta_pred <- sqrt(exp(pred_mean[, 2]))
# sigma_eps_pred <- sqrt(exp(pred_mean[, 3]))

# # test_sample <- test_output[sample, ] # arctanh(phi), log(sigma_eta), log(sigma_eps)
# phi_true <- tanh(test_output[, 1])
# sigma_eta_true <- sqrt(exp(test_output[, 2]))
# sigma_eps_true <- sqrt(exp(test_output[, 3]))

# phi_df <- data.frame(param = "phi", pred = phi_pred, true = phi_true)
# sigma_eta_df <- data.frame(param = "sigma[eta]", pred = sigma_eta_pred, true = sigma_eta_true)
# sigma_eps_df <- data.frame(param = "sigma[epsilon]", pred = sigma_eps_pred, true = sigma_eps_true)
# mean_pred_df <- rbind(phi_df, sigma_eta_df, sigma_eps_df)

# mean_pred_plot <- ggplot(mean_pred_df, aes(x = true, y = pred)) + 
#   geom_point(color = "salmon") + 
#   geom_abline(intercept = 0, slope = 1, color = "black") +
#   facet_wrap(~param, scales = "free", labeller = label_parsed) +
#   theme(aspect.ratio = 1) +
#   # guides(fill="none") +
#   xlab("True parameters") +
#   ylab("Predicted parameters") +
#   # coord_fixed() +
#   # scale_x_continuous(expand = c(0, 0)) + 
#   # scale_y_continuous(expand = c(0, 0)) +
#   theme_bw() +
#   theme(text = element_text(size = 30))

# png(paste0(output_dir, "/lgss_pred_vs_test_plot_", data_date, ".png"), width = 1200, height = 400)
# print(mean_pred_plot)  
# dev.off()

# ## Variance of the predicted parameters
# construct_L_mat <- function(v, d) {
#   L_mat <- diag(exp(v[1:d]))
#   diag(L_mat[-1, -d]) <- v[(d + 1):length(v)]
#   return(L_mat)
# }

# pred_chol <- pred_output[, (n_mean_elements + 1):output_dim]
# pred_chol_ls <- lapply(1:nrow(pred_chol), function(i) construct_L_mat(pred_chol[i, ], n_mean_elements))
# L_inv_ls <- lapply(pred_chol_ls, solve)

## For each test sample, simulate from the predicted posterior distribution
S <- 10000

phi_mean <- c()
phi_quantiles <- matrix(NA, nrow = length(test_output), ncol = 2)
phi_post_samples_ls <- list()
for (s in 1:length(test_output)) {
    # post_samples <- rmvnorm(S, pred_mean[s, ], t(L_inv_ls[[s]]) %*% L_inv_ls[[s]])
    post_samples <- rnorm(S, pred_mean[s], pred_sd[s])
    
    ### Back-transform
    if (use_arctanh) {
        phi_post_samples <- tanh(post_samples)
    } else { # use logit
        phi_post_samples <- exp(post_samples) / (1 + exp(-post_samples))
    }
    phi_mean[s] <- mean(phi_post_samples)
    phi_quantiles[s, ] <- quantile(phi_post_samples, c(0.025, 0.975))
    phi_post_samples_ls[[s]] <- phi_post_samples
}

if (use_arctanh) {
    phi_true <- tanh(test_output)
} else { # use logit
    phi_true <- exp(test_output) / (1 + exp(-test_output))
}

## Then maybe plot all of the samples point-wise
phi_df <- data.frame(pred = phi_mean, true = phi_true,
                     lower = phi_quantiles[, 1], upper = phi_quantiles[, 2])

phi_coverage <- sum(phi_df$lower <= phi_df$true & phi_df$true <= phi_df$upper) / nrow(phi_df)

## Plot the mean against true values
mean_plot <- phi_df %>%
              ggplot() +
              geom_point(aes(x = true, y = pred), color = "salmon") +
              geom_abline(intercept = 0, slope = 1, color = "black") +
              # xlab(expression("True "~phi)) +
              # ylab(expression("Posterior mean of "~phi)) +
              xlab("True values") +
              ylab("Posterior means") +
              theme_bw() +
              theme(text = element_text(size = 24))

## Plot the mean against the true values
png(paste0("./plots/lgss_mean_plot_", data_date, ".png"), width = 500, height = 500)
# plot(phi_df$true, phi_df$pred, xlab = "True phi", ylab = "Predicted phi", col = "salmon")
# abline(a = 0, b = 1, col = "black")
print(mean_plot)
dev.off()

set.seed(2025)
n_plot_samples <- 100
phi_quantile_plot <- phi_df %>% 
    slice_sample(n = n_plot_samples) %>%
    ggplot() + 
    geom_errorbar(aes(x = 1:n_plot_samples, ymin = lower, ymax = upper), width = 0, color = "black", lwd = 1) +
    geom_point(aes(x = 1:n_plot_samples, y = pred), color = "black", size = 3, alpha = 0.5) + 
    geom_point(aes(x = 1:n_plot_samples, y = true), color = "salmon", size = 3, alpha = 0.75) +
    xlab("Test sample") +
    ylab(bquote(phi)) +
    theme_bw() +
    theme(text = element_text(size = 24)) 


png(paste0("./plots/lgss_quantile_plot_", data_date, ".png"), width = 1000, height = 300)
print(phi_quantile_plot)
dev.off()


## Compare to output from HMC
set.seed(2025)
test_sample <- sample(1:length(test_output), 5)

hmc_vs_cnn_plots <- list()

for (i in 1:length(test_sample)) {
  s <- test_sample[i]
  hmc_output <- qread(paste0("output/lgss_hmc_results_", 
                            formatC(s, width = 2, format = "d", flag = "0"), 
                            "_", data_date, ".qs"))
  hmc.phi <- as.vector(hmc_output$draws[-(1:hmc_output$burn_in), , 1])

  hmc_vs_cnn_df <- data.frame(hmc = hmc.phi, cnn = phi_post_samples_ls[[s]])
  phi_true_df <- data.frame(phi_true = phi_true[s])

  ## Plot HMC vs CNN posterior
  hmc_vs_cnn_plot <- ggplot(hmc_vs_cnn_df) +
    geom_density(aes(x = hmc, col = "HMC"), lwd = 1) +
    geom_density(aes(x = cnn, col = "CNN"), lwd = 1) +
    geom_vline(data = phi_true_df, aes(xintercept = phi_true), lty = 2, lwd = 1) +
    xlab(bquote(phi)) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    guides(color = "none")
    # labs(col = "Method")

  hmc_vs_cnn_plots[[i]] <- hmc_vs_cnn_plot

  
}


png(paste0("./plots/lgss_hmc_vs_cnn_", data_date, ".png"), width = 1500, height = 400)
grid.arrange(grobs = hmc_vs_cnn_plots, ncol = length(test_sample))
dev.off()




