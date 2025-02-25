## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/phi/")

rm(list = ls())

library(mvtnorm)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)
library(qs)

source("./source/create_model.R")
source("./source/posterior_loss.R")


# #List physical devices
# gpus <- tf$config$experimental$list_physical_devices('GPU')

# if (length(gpus) > 0) {
#   tryCatch({
#     # Restrict TensorFlow to only allocate 4GB of memory on the first GPU
#     tf$config$experimental$set_virtual_device_configuration(
#       gpus[[1]],
#       list(tf$config$experimental$VirtualDeviceConfiguration(memory_limit=4096*10))
#     )
    
#     logical_gpus <- tf$config$experimental$list_logical_devices('GPU')
    
#     print(paste0(length(gpus), " Physical GPUs,", length(logical_gpus), " Logical GPUs"))
#   }, error = function(e) {
#     # Virtual devices must be set before GPUs have been initialized
#     print(e)
#   })
# }

## Read data
data_date <- "20250108"
simulated_data <- qread(paste0("./data/lgss_data_", data_date, ".qs"))
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
if (!dir.exists(output_dir)) {
  dir.create(paste0(output_dir))
} #else { # delete all previously saved checkpoints
#   unlink(paste0(output_dir, "/*"))
# }

checkpoint_path <- paste0(output_dir, "checkpoints/cp-{epoch:04d}.ckpt")

batch_size <- 64
epochs <- 20

## Train the model
print("Training CNN...")
# checkpoint_path <- paste0("output/", setsvf, "/checkpoints/cp-{epoch:04d}.ckpt")
# checkpoint_dir <- fs::path_dir(checkpoint_path)

cp_callback <- callback_model_checkpoint(
  filepath = checkpoint_path,
  save_weights_only = TRUE,
  verbose = 1#,m
  # save_freq = 10*batch_size # save every 10 epochs
)

# Train the model with the new callback
history <- model %>% fit(
    train_input, 
    train_output,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(val_input, val_output),
    callbacks = list(cp_callback)
)

qsave(history, file = paste0(output_dir, "history_", data_date, ".qs"))

png(paste0("./plots/loss_plot_", data_date, ".png"))
plot(history)
dev.off()

# ## Choose model with lowest validation loss
# checkpt <- which.min(history$metrics$val_loss) 
# cp_restart <- paste0(output_dir, "checkpoints/cp-", formatC(checkpt, width = 4, flag = "0"), ".ckpt")
# # latest <- tf$train$latest_checkpoint(checkpoint_dir)
# load_model_weights_tf(model, cp_restart)

# ## Predict on test set
# pred_time <- system.time({
#     pred_output <- model %>% predict(test_input)
# })
# cat("Prediction time: ", pred_time[3], "\n")

# pred_mean <- pred_output[, 1:n_mean_elements]

# # Compare predicted mean vs true parameters for a chosen sample in the test set
# phi_pred <- tanh(pred_mean[, 1]) # back-transform
# sigma_eta_pred <- sqrt(exp(pred_mean[, 2]))
# sigma_eps_pred <- sqrt(exp(pred_mean[, 3]))

# # test_sample <- test_output[sample, ] # arctanh(phi), log(sigma_eta), log(sigma_eps)
# phi_true <- tanh(test_output[, 1])
# sigma_eta_true <- sqrt(exp(test_output[, 2]))
# sigma_eps_true <- sqrt(exp(test_output[, 3]))

# phi_df <- data.frame(param = "phi", pred = phi_pred, true = phi_true)
# sigma_eta_df <- data.frame(param = "sigma_eta", pred = sigma_eta_pred, true = sigma_eta_true)
# sigma_eps_df <- data.frame(param = "sigma_eps", pred = sigma_eps_pred, true = sigma_eps_true)
# mean_pred_df <- rbind(phi_df, sigma_eta_df, sigma_eps_df)

# png(paste0(output_dir, "/pred_vs_test_plot_", data_date, ".png"), width = 1200, height = 400)

# ggplot(mean_pred_df, aes(x = true, y = pred, color = param)) + 
#   geom_point() + 
#   geom_abline(intercept = 0, slope = 1, color = "black") +
#   facet_wrap(~param, scales = "free") +
#   theme(aspect.ratio = 1) +
#   theme_bw() +
#   theme(text = element_text(size = 20)) 
  
# dev.off()

# ## Variance of the predicted parameters
# construct_L_mat <- function(v, d) {
#   L_mat <- diag(exp(v[1:d]))
#   diag(L_mat[-1, -d]) <- v[(d + 1):length(v)]
#   return(L_mat)
# }

# pred_chol <- pred_output[, (n_mean_elements + 1):output_dim]
# pred_chol_ls <- lapply(1:nrow(pred_chol), function(i) construct_L_mat(pred_chol[i, ], n_mean_elements))

# ## For each test sample, simulate from the predicted posterior distribution
# s <- 1
# S <- 1000
# post_samples <- rmvnorm(S, pred_mean[s, ], pred_chol_ls[[s]] %*% t(pred_chol_ls[[s]]))

# ### Back-transform
# phi_post_samples <- tanh(post_samples[, 1])
# sigma_eta_post_samples <- sqrt(exp(post_samples[, 2]))
# sigma_eps_post_samples <- sqrt(exp(post_samples[, 3]))

# ## Need to obtain the quantiles here

# ## Then maybe plot all of the samples point-wise

# png(paste0(output_dir, "/post_samples_plot_", data_date, ".png"))

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
