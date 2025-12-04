## Train a CNN on ice thickness and velocity data

setwd("~/SSA_model/CNN/real_data/")

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
# output_var <- "all" # "all" #"bed"  # "grounding_line" # "bed"
use_missing_pattern <- T
leave_one_out <- T

source("./source/create_model.R")
# source("./source/custom_loss_function.R")
source("./source/posterior_loss.R")

## Read data
data_date <- "20241111" #"20241103" #"20241103"
# arg <- commandyArgs(trailingOnly = TRUE)
sets <- 51:100 #c(1,3,5) #11:15 #6:10 #arg
# setf <- formatC(set, width=2, flag="0")
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Directories
train_data_dir <- paste0("./data/training_data", "/", setsf, "/")
output_dir <- paste0("./output/cnn/", setsf, "/")
plot_dir <- paste0("./plots/cnn/", setsf, "/")

if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
}

print("Reading data...")
system.time({
  train_data <- qread(file = paste0(train_data_dir, "train_data_", data_date, ".qs"))
  val_data <- qread(file = paste0(train_data_dir, "val_data_", data_date, ".qs"))
  test_data <- qread(file = paste0(train_data_dir, "test_data_", data_date, ".qs"))
})

train_input <- train_data$input
val_input <- val_data$input
test_input <- test_data$input

if (leave_one_out) {
    n_years <- dim(train_input)[3] - 1
    train_input <- train_input[,,1:n_years,]
    val_input <- val_input[,,1:n_years,]
    test_input <- test_input[,,1:n_years,]
    print(paste0("Using first ", n_years, " years of data for leave-one-out cross-validation"))
}

train_output <- cbind(train_data$fric_coefs, train_data$bed_coefs, train_data$grounding_line)
val_output <- cbind(val_data$fric_coefs, val_data$bed_coefs, val_data$grounding_line)
test_output <- cbind(test_data$fric_coefs, test_data$bed_coefs, test_data$grounding_line)

png(paste0(plot_dir, "input_", data_date, ".png"), width = 800, height = 400)
par(mfrow = c(1, 2))
sp <- sample(1:dim(val_input)[1], 1)
# matplot(train_input[100,,,1], type = "l", col = "grey")
# matlines(val_input[100,,,1], col = scales::alpha("red", 0.1))
image(train_input[sp,, ,1], col = terrain.colors(100), main = "Train input example")
image(val_input[sp,, ,1], col = terrain.colors(100), main = "Val input example")
dev.off()

sp <- sample(1:dim(val_output)[1], 1)
png(paste0(plot_dir, "output_", data_date, ".png"))
plot(train_output[sp, ], type = "l", col = "grey", 
    main = "Output example", ylab = "Output value", xlab = "Coefficient index")
lines(val_output[sp, ], col = scales::alpha("red", 0.5))
legend("topright", legend = c("train", "val"), col = c("grey", "red"), lty = 1)
dev.off()

############################
##      Create model      ##
############################
print("Creating model...")
# Create a basic model instance
input_dim <- dim(train_input)[2:4]
n_fric_basis <- dim(train_data$fric_coefs)[2]
n_bed_basis <- dim(train_data$bed_coefs)[2]
n_gl <- dim(train_data$grounding_line)[2]
n_mean_elements <- n_fric_basis + n_bed_basis + n_gl
n_chol_elements <- n_mean_elements + (n_mean_elements - 3) # diagonal + lower-diag elements
output_dim <- n_mean_elements + n_chol_elements  # THIS NEEDS TO CHANGE TO THE TOTAL NUMBER OF BASIS FUNCTIONS + COVARIANCE PARAMETERS

  model <- create_model_posterior(input_dim = input_dim, 
                                  output_dim = output_dim,
                                  n_bed_basis = n_bed_basis,
                                  n_fric_basis = n_fric_basis,
                                  n_gl = n_gl)

# Display the model's architecture
summary(model)

# Compute training loss (and metrics) before training
train_metrics <- model %>% evaluate(train_input, train_output, verbose = 0)

# Compute validation loss (and metrics) before training
val_metrics <- model %>% evaluate(val_input, val_output, verbose = 0)

cat("Initial training loss:", train_metrics["loss"], "\n")
cat("Initial validation loss:", val_metrics["loss"], "\n")

if (!dir.exists(output_dir)) {
  dir.create(paste0(output_dir))
} else { # delete all previously saved checkpoints
  unlink(paste0(output_dir, "/*"))
}

checkpoint_path <- paste0(output_dir, "checkpoints/cp-{epoch:04d}.ckpt")
# checkpoint_dir <- fs::path_dir(checkpoint_path)

batch_size <- 64
epochs <- 50

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

  
  qsave(history, file = paste0(output_dir, "/history_", data_date, ".qs"))

  # Save the entire model as a SavedModel.
  save_model_tf(model, paste0(output_dir, "/model_", data_date))
} else {
  model <- load_model_tf(paste0(output_dir, "/model_", data_date))
  history <- qread(file = paste0(output_dir, "/history_", data_date, ".qs"))
}

## Plot the loss
loss_plot <- history %>%
  plot() +
  coord_cartesian(xlim = c(1, epochs))

png(paste0(plot_dir, "loss_", data_date, ".png"))
print(loss_plot)
dev.off()
