## data generator for CNN

setwd("/home/babv971/SSA_model/CNN/pilot/")
rm(list = ls())

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)

source("./source/create_model.R")
# generate_data <- function(input_path, output_path) {
#   training_data <- data_generator(input_path, output_path)

#   return(training_data)
# }


# data_generator <- function(input_paths, output_paths) {

#   index <- 1
  
#   function() {
#     while(TRUE) {
#       if (index > length(file_paths)) {
#         index <<- 1
#       }

#         thickness_velocity <- readRDS(input_paths[[index]])
#         friction <- readRDS(output_paths[[index]])

#         input <- thickness_velocity
#         output <- friction$basis_coefs

#         index <<- index + 1
#         list(input, output)

#     }
#   }
# }



# Define the generator functions
read_rds_batch <- function(file_paths, batch_size) {
  index <- 1
  
  function() {
    while(TRUE) {
      if (index > length(file_paths)) {
        index <<- 1
      }
      
      train_data <- readRDS(file_paths[index])
      index <<- index + 1
      
      features <- train_data$input #as.matrix(data[, -ncol(data)])
      labels <- train_data$output #as.matrix(data[, ncol(data)])
      
      # for (start in seq(1, nrow(features), by = batch_size)) {
      #   end <- min(start + batch_size - 1, nrow(features))
      #   yield(list(features[start:end, , drop = FALSE], labels[start:end, , drop = FALSE]))
      # }

      list(features, labels)
    }
  }
}

data_generator <- function(file_paths, batch_size) {
  generator <- read_rds_batch(file_paths, batch_size)
  
  function() {
    batch <- generator()
    list(as_tensor(batch[[1]]), as_tensor(batch[[2]])) # why use so many functions???
  }
}

# Prepare file paths and batch size
# file_paths <- list.files("path/to/your/rds/files", pattern = "\\.rds$", full.names = TRUE)

data_date <- "20240320" #"20220329" 
sets <- 1:2
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Work out how many simulations there are
n_per_set <- 1000
n_sims <- n_per_set * length(sets)

## Work out how many observations are needed for a 89/10/1 split
n_train <- floor(n_sims * 0.89)

last_train_set <- ceiling(n_train / n_per_set)

train_sets <- 1:last_train_set
file_paths <- lapply(train_sets, function(x) paste0("./training_data/", setsf, "/train_data_", formatC(x, width=2, flag="0"), "_", data_date, ".rds"))

batch_size <- 64

# Calculate steps per epoch
steps_per_epoch <- n_train/batch_size

# test <- generate_data(input_paths, output_paths)

# # Define and compile your model
# n_features <- 10  # Adjust to match your data

# model <- keras_model_sequential() %>%
#   layer_dense(units = 128, activation = 'relu', input_shape = c(n_features)) %>%
#   layer_dense(units = 64, activation = 'relu') %>%
#   layer_dense(units = 1, activation = 'sigmoid')

# model %>% compile(
#   optimizer = 'adam',
#   loss = 'binary_crossentropy',
#   metrics = c('accuracy')
# )

output_dim <- 150 #ncol(train_output)
model <- create_model(output_dim = output_dim)

# Create the data generator and train the model
train_generator <- data_generator(file_paths, batch_size)

model %>% fit(
  train_generator,
  steps_per_epoch = steps_per_epoch,
  epochs = 10
)
