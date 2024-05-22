## Train a CNN on ice thickness and velocity data

setwd("/home/babv971/SSA_model/CNN/pilot/")

rm(list = ls())

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)

source("./source/create_model.R")
source("./source/data_generator.R")

# List physical devices
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

## Flags
# rerun_cnn <- T
# save_output <- T

## Read data
data_date <- "20240320"
# arg <- commandArgs(trailingOnly = TRUE)
sets <- 1:10 #arg
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

print("Reading data...")
train_data <- readRDS(file = paste0("./training_data/", setsf, "/train_data_", data_date, ".rds"))
val_data <- readRDS(file = paste0("./training_data/", setsf, "/val_data_", data_date, ".rds"))
test_data <- readRDS(file = paste0("./training_data/", setsf, "/test_data_", data_date, ".rds"))

train_input <- train_data$input
val_input <- val_data$input
test_input <- test_data$input

train_output <- train_data$output
val_output <- val_data$output
test_output <- test_data$output

# if (rerun_cnn) {
print("Training CNN...")
# Create a basic model instance
output_dim <- ncol(train_output)
model <- create_model(output_dim = output_dim)

# Display the model's architecture
summary(model)

# Create a callback that saves the model's weights
output_dir <- paste0("./output/", setsf) #, "/checkpoints")

if (file.exists(output_dir)) { # remove old checkpoints
    unlink(output_dir, recursive = TRUE)
  } else {
    dir.create(output_dir)
  }

checkpoint_path <- paste0("output/", setsf, "/checkpoints/cp-{epoch:04d}.ckpt")
# checkpoint_dir <- fs::path_dir(checkpoint_path)

batch_size <- 64
epochs <- 10

input_paths <- lapply(sets, function(x) paste0("./training_data/thickness_velocity_arr_", 
                                                formatC(x, width=2, flag="0"), "_", data_date, ".rds"))
output_paths <- lapply(sets, function(x) paste0("./training_data/friction_basis_", 
                                                formatC(x, width=2, flag="0"), "_", data_date, ".rds"))

# if (rerun_cnn) {

  cp_callback <- callback_model_checkpoint(
    filepath = checkpoint_path,
    save_weights_only = TRUE,
    verbose = 1#,
    # save_freq = 10*batch_size # save every 10 epochs
  )

  # Train the model with the new callback
  history <- model %>% fit(
      # train_input, 
      # train_output,
      data_generator(input_paths, output_paths),
      epochs = epochs,
      batch_size = batch_size,
      validation_data = list(val_input, val_output),
      callbacks = list(cp_callback)
  )

  
  saveRDS(history, file = paste0("./output/", setsf, "/history_", data_date, ".rds"))

  # Save the entire model as a SavedModel.
  save_model_tf(model, paste0("output/", setsf, "/model_", data_date))
# } else {
#   model <- load_model_tf(paste0("output/", setsf, "/model_", data_date))
#   history <- readRDS(file = paste0("./output/", setsf, "/history_", data_date, ".rds"))
# }

