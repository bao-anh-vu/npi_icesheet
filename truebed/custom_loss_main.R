## Example from https://towardsdatascience.com/custom-loss-functions-for-deep-learning-predicting-home-values-with-keras-for-r-532c9e098d1f

source("~/SSA_model/CNN/truebed/source/custom_loss_function.R")

# load the data set
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)

data <- dataset_boston_housing()
c(c(train_data,train_targets), c(test_data, test_targets)) %<-% data

# transform the training and test labels
train_targets <- (train_targets*1000)^2/2500000
test_targets <- (test_targets*1000)^2/2500000

# The model as specified in "Deep Learning with R"
model <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = "relu",
              input_shape = dim(train_data)[[2]]) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 1)

# Compile the model, and select one of the loss functions
q <- 0.95 #c(0.1, 0.5, 0.9)

# losses <- c(keras::loss_mean_squared_error,  
#     keras::loss_mean_squared_logarithmic_error, 
#     MLAE, MSLAE, quantile_loss_wrap)

# metrics <- c("mae", "mse", "msle", "mlae", "mslae")    
model %>% compile(
  optimizer = "rmsprop",
  loss = quantile_loss_wrap(q), #losses[1],
  metrics = quantile_loss_wrap(q) #c("mae") # is this the right thing to do?
)

# Train the model with validation
model %>% fit(
  train_data,
  train_targets,
  epochs = 20,
  batch_size = 5,
  verbose = 1,
  validation_split = 0.2
)

# Calculate the mean absolute error 
results <- model %>% evaluate(test_data, test_targets, verbose = 0)
results

