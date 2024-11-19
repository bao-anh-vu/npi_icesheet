## Train a CNN on ice thickness and velocity data

setwd("~/SSA_model/CNN/simbed/")

rm(list = ls())

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

## Flags
rerun_cnn <- T
sim_beds <- T
output_var <- "all" # "all" #"bed"  # "grounding_line" # "bed"
use_missing_pattern <- F
# save_output <- T

source("./source/create_model.R")
# source("./source/custom_loss_function.R")
source("./source/posterior_loss.R")

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
sets <- 1:50 #c(1,3,5) #11:15 #6:10 #arg
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

print("Reading data...")
if (use_missing_pattern) {
  train_data_dir <- paste0("./training_data", "/", setsf, "/missing")
} else {
  train_data_dir <- paste0("./training_data", "/", setsf, "/nonmissing")
}
system.time({
  train_data <- qread(file = paste0(train_data_dir, "/train_data_", data_date, ".qs"))
  val_data <- qread(file = paste0(train_data_dir, "/val_data_", data_date, ".qs"))
  test_data <- qread(file = paste0(train_data_dir, "/test_data_", data_date, ".qs"))
})
# system.time({
  # qsave(train_data, file = paste0(train_data_dir, "/train_data_", data_date, ".qs"))
  # qsave(val_data, file = paste0(train_data_dir, "/val_data_", data_date, ".qs"))
  # qsave(test_data, file = paste0(train_data_dir, "/test_data_", data_date, ".qs"))
# })

train_input <- train_data$input
val_input <- val_data$input
test_input <- test_data$input

train_output <- cbind(train_data$fric_coefs, train_data$bed_coefs, train_data$grounding_line)
val_output <- cbind(val_data$fric_coefs, val_data$bed_coefs, val_data$grounding_line)
test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)


# train_output <- train_data$output
# val_output <- val_data$output
# test_output <- test_data$output

# if (rerun_cnn) {
print("Creating model...")
# Create a basic model instance
input_dim <- dim(train_data$input)[2:4]
n_basis_funs <- dim(train_data$fric_coefs)[2]
n_gl <- dim(train_data$grounding_line)[2]
n_mean_elements <- n_basis_funs * 2 + n_gl
n_chol_elements <- (n_basis_funs * 2 + n_gl) + (n_basis_funs - 1) * 2 + (n_gl - 1) # diagonal + lower-diag elements
output_dim <- n_mean_elements + n_chol_elements  # THIS NEEDS TO CHANGE TO THE TOTAL NUMBER OF BASIS FUNCTIONS + COVARIANCE PARAMETERS


# if (output_var == "friction") {
#   model <- create_model(input_dim = input_dim, output_dim = output_dim)
# } else if (output_var == "grounding_line") {
#   model <- create_model(input_dim = input_dim, output_dim = output_dim)
# } else if (output_var == "bed") { ## bed
#   model <- create_model_bed(input_dim = input_dim, output_dim = output_dim)
# } else { ## all variables
  model <- create_model_posterior(input_dim = input_dim, 
                                  output_dim = output_dim,
                                  n_basis_funs = n_basis_funs,
                                  n_gl = n_gl)
# }

# Display the model's architecture
summary(model)

# Create a callback that saves the model's weights
if (use_missing_pattern) {
  output_dir <- paste0("./output/posterior/", setsf, "/missing")
} else {
  output_dir <- paste0("./output/posterior/", setsf, "/nonmissing")
}

if (!dir.exists(output_dir)) {
  dir.create(paste0(output_dir))
} else { # delete all previously saved checkpoints
  unlink(paste0(output_dir, "/*"))
}

checkpoint_path <- paste0(output_dir, "/checkpoints/cp-{epoch:04d}.ckpt")
# checkpoint_dir <- fs::path_dir(checkpoint_path)

batch_size <- 64
epochs <- 20

if (rerun_cnn) {
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

  
  saveRDS(history, file = paste0(output_dir, "/history_", data_date, ".rds"))

  # Save the entire model as a SavedModel.
  save_model_tf(model, paste0(output_dir, "/model_", data_date))
} else {
  model <- load_model_tf(paste0(output_dir, "/model_", data_date))
  history <- readRDS(file = paste0(output_dir, "/history_", data_date, ".rds"))
}

loss_plot <- history %>%
  plot() +
  coord_cartesian(xlim = c(1, epochs))
# ## Plot the loss
png(paste0("plots/posterior/loss_", data_date, ".png"))
print(loss_plot)
dev.off()
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

