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
# n_chol_elements <- 5 # diagonal + lower-diag elements
output_dim <- n_mean_elements #+ n_chol_elements  

model <- create_model(input_dim = input_dim, 
                      output_dim = output_dim) 

# Display the model's architecture
summary(model)

output_dir <- "./output_sq_error/"
if (!dir.exists(output_dir)) {
  dir.create(paste0(output_dir))
} #else { # delete all previously saved checkpoints
#   unlink(paste0(output_dir, "/*"))
# }

checkpoint_path <- paste0(output_dir, "/checkpoints/cp-{epoch:04d}.ckpt")

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

qsave(history, file = paste0(output_dir, "/history_", data_date, ".qs"))

png(paste0(output_dir, "/loss_plot_", data_date, ".png"))
plot(history)
dev.off()

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

png(paste0(output_dir, "/pred_vs_test_plot_", data_date, ".png"))
par(mfrow = c(3, 1))
plot(phi_pred, phi_true)
abline(0, 1, col = "red")

plot(sigma_eta_pred, sigma_eta_true)
abline(0, 1, col = "red")

plot(sigma_eps_pred, sigma_eps_true)
abline(0, 1, col = "red")
dev.off()


# ## Predict on the ``real'' data
# y <- qread(paste0("output/y_", data_date, ".qs"))
# y_std <- (y - mean_input) / sd_input
# y_std <- array_reshape(y_std, c(1, length(y_std), 1))

# real_pred <- model %>% predict(y_std)
# real_pred_mean <- real_pred[, 1:n_mean_elements]

# phi_pred <- tanh(real_pred_mean[1])
# sigma_eta_pred <- exp(real_pred_mean[2])
# sigma_eps_pred <- exp(real_pred_mean[3])
