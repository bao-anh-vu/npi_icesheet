## Train a CNN on ice thickness and velocity data

setwd("~/SSA_model/CNN/pilot/simbed/")

rm(list = ls())

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ggplot2)

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
rerun_cnn <- T
sim_beds <- T
output_var <- "all" #"bed" # "friction" # "grounding_line" # "bed"
# save_output <- T

source("./source/create_model.R")

# if (output_var == "friction") {
#   source("./source/create_model.R")
# } else if (output_var == "grounding_line") {
#   source("./source/create_cnn_gl.R")
# } else if (output_var == "bed") {
#   source("./source/create_model.R")
# }

## Read data
data_date <- "20240320"
# arg <- commandArgs(trailingOnly = TRUE)
sets <- 1:2 #arg
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

print("Reading data...")
# if (sim_beds) {
  train_data_dir <- "~/SSA_model/CNN/pilot/training_data_bed"
  train_data <- readRDS(file = paste0(train_data_dir, "/", setsf, "/train_data_", data_date, ".rds"))
  val_data <- readRDS(file = paste0(train_data_dir, "/", setsf, "/val_data_", data_date, ".rds"))
  test_data <- readRDS(file = paste0(train_data_dir, "/", setsf, "/test_data_", data_date, ".rds"))

# } else {
#   train_data_dir <- "./training_data"
#   train_data <- readRDS(file = paste0(train_data_dir, "/", output_var, "/", setsf, "/train_data_", data_date, ".rds"))
#   val_data <- readRDS(file = paste0(train_data_dir, "/", output_var, "/", setsf, "/val_data_", data_date, ".rds"))
#   test_data <- readRDS(file = paste0(train_data_dir, "/", output_var, "/", setsf, "/test_data_", data_date, ".rds"))

# }

train_input <- train_data$input
val_input <- val_data$input
test_input <- test_data$input

if (output_var == "friction") {
  train_output <- train_data$fric_coefs
  val_output <- val_data$fric_coefs
  test_output <- test_data$fric_coefs
} else if (output_var == "grounding_line") {
  train_output <- train_data$grounding_line
  val_output <- val_data$grounding_line
  test_output <- test_data$grounding_line
} else if (output_var == "bed") {
  train_output <- train_data$bed_coefs
  val_output <- val_data$bed_coefs
  test_output <- test_data$bed_coefs
} else if (output_var == "all") {
  train_output <- cbind(train_data$fric_coefs, train_data$bed_coefs)
  val_output <- cbind(val_data$fric_coefs, val_data$bed_coefs)
  test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs)
}

# train_output <- train_data$output
# val_output <- val_data$output
# test_output <- test_data$output

# if (rerun_cnn) {
print("Training CNN...")
# Create a basic model instance
output_dim <- ncol(train_output)

if (output_var == "friction") {
  model <- create_model(output_dim = output_dim)
} else if (output_var == "grounding_line") {
  model <- create_model(output_dim = output_dim)
} else if (output_var == "bed") { ## bed
  model <- create_model_bed(output_dim = output_dim)
} else {
  model <- create_model(output_dim = output_dim)
}

# Display the model's architecture
summary(model)

# Create a callback that saves the model's weights
# if (sim_beds) {
  output_dir <- paste0("./output/", output_var, "/", setsf)
# } else {
#   output_dir <- paste0("./output/", output_var, "/", setsf)
# }

if (!dir.exists(output_dir)) {
  dir.create(paste0(output_dir))
} else { # delete all previously saved checkpoints
  unlink(paste0(output_dir, "/*"))
}
# dir.create(paste0("output/", setsf)
checkpoint_path <- paste0(output_dir, "/checkpoints/cp-{epoch:04d}.ckpt")
# checkpoint_dir <- fs::path_dir(checkpoint_path)

batch_size <- 64
epochs <- 10

if (rerun_cnn) {
  
  # checkpoint_path <- paste0("output/", setsf, "/checkpoints/cp-{epoch:04d}.ckpt")
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

  
  saveRDS(history, file = paste0(output_dir, "/history_", data_date, ".rds"))

  # Save the entire model as a SavedModel.
  save_model_tf(model, paste0(output_dir, "/model_", data_date))
} else {
  model <- load_model_tf(paste0(output_dir, "/model_", data_date))
  history <- readRDS(file = paste0(output_dir, "/history_", data_date, ".rds"))
}

# ## Plot the loss

history %>%
  plot() +
  coord_cartesian(xlim = c(1, epochs))

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

