## Train a CNN on ice thickness and velocity data

setwd("/home/babv971/SSA_model/CNN/")

rm(list = ls())

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)

source("source/create_model.R")
## Flags
rerun_cnn <- T
# save_output <- T

## Read data
data_date <- "20240320"
# arg <- commandArgs(trailingOnly = TRUE)
sets <- 1 #arg
# setf <- formatC(set, width=2, flag="0")
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

print("Reading data...")
train_data <- readRDS(file = paste0("./training_data/train_data_", setf, "_", data_date, ".rds"))
train_input <- train_data$train_input
val_input <- train_data$val_input
test_input <- train_data$test_input

train_output <- train_data$train_output
val_output <- train_data$val_output
test_output <- train_data$test_output

# if (rerun_cnn) {
print("Training CNN...")
# Create a basic model instance
model <- create_model()

# Display the model's architecture
summary(model)

# Create a callback that saves the model's weights
checkpoint_path <- paste0("output/checkpoints", setf, "/:/cp-{epoch:04d}.ckpt")
checkpoint_dir <- fs::path_dir(checkpoint_path)

batch_size <- 64
epochs <- 100

if (rerun_cnn) {
  cp_callback <- callback_model_checkpoint(
    filepath = checkpoint_path,
    save_weights_only = TRUE,
    verbose = 1#,
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

  
  saveRDS(history, file = paste0("./output/history_", setf, "_", data_date, ".rds"))

  # Save the entire model as a SavedModel.
  save_model_tf(model, paste0("output/model_", setf, "_", data_date))
} else {
  model <- load_model_tf(paste0("output/model_", setf, "_", data_date))
  history <- readRDS(file = paste0("./output/history_", setf, "_", data_date, ".rds"))
}

# ## Plot the loss

# history %>%
#   plot() +
#   coord_cartesian(xlim = c(1, epochs))

# ## Get rid of first training loss
# plot(history$metrics$loss[2:60],type = "l")
# lines(history$metrics$val_loss[2:60], col = "red")

# browser()

# ## Predict friction coefficients on test set
# pred_coefs_new <- model %>% predict(test_input)
# saveRDS(pred_coefs_new, file = paste0("./output/pred_coefs_", setf, "_", data_date, ".rds"))

# results <- model %>% evaluate(test_input, test_output, batch_size = batch_size)
# cat("test loss, test acc:", results)

# # Save the entire model as a SavedModel.
# save_model_tf(model, "output/my_model")

# restored_model <- load_model_tf('output/my_model')

# # Re-evaluate the model
# new <- restored_model %>% fit(
#     train_input, 
#     train_output,
#     epochs = 10,
#     batch_size = 64,
#     validation_data = list(val_input, val_output),
#     callbacks = list(cp_callback)
# )

