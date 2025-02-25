## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/phi")

source("./source/sim_data.R")

library(mvtnorm)
library(tensorflow)
library(qs)

# #List physical devices
# gpus <- tf$config$experimental$list_physical_devices('GPU')

# if (length(gpus) > 0) {
#   tryCatch({
#     # Restrict TensorFlofw to only allocate 4GB of memory on the first GPU
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
regenerate_data <- T
save_data <- T
use_arctanh <- T # use either arctanh or logit

## Generate data
data_date <- "20250108"

# sigma_eps <- 0.5 # measurement error var
# sigma_eta <- 0.7 # process error var
# phi <- 0.9

# Generate "real" process x_{1:T} and data y_{1:T}
iters <- 1000 #T

## Plot some examples to make sure the data is being generated correctly
test1 <- sim_data(phi = 0.99, iters = iters)
test2 <- sim_data(phi = 0.9, iters = iters)
test3 <- sim_data(phi = 0.1, iters = iters)
test4 <- sim_data(phi = -0.9, iters = iters)

# test1 <- sim_data(x_ini = x_ini, phi = 0.9, sigma_eta = 2.0, sigma_eps = 0.5, iters = iters)
# test2 <- sim_data(x_ini = x_ini, phi = 0.9, sigma_eta = 1.0, sigma_eps = 0.5, iters = iters)
# test3 <- sim_data(x_ini = x_ini, phi = 0.9, sigma_eta = 0.5, sigma_eps = 0.5, iters = iters)
# test4 <- sim_data(x_ini = x_ini, phi = 0.9, sigma_eta = 0.1, sigma_eps = 0.5, iters = iters)


plot_range <- 1:100
png("./plots/lgss_compare_y.png", width = 800, height = 500)
plot(test1$y[plot_range], type = "l")
lines(test2$y[plot_range], col = "red")
lines(test3$y[plot_range], col = "salmon")
lines(test4$y[plot_range], col = "blue")
dev.off()


## Simulate training data
nsims <- 50000#0
prior_mean <- 0
prior_sd <- 1

if (regenerate_data) {
  print("Simulating data...")
  ## Prior for the parameters: arctanh(phi) ~ N(0, 1), log(sigma_eta) ~ N(0, 1), log(sigma_eps) ~ N(0, 1)
  param_samples <- rnorm(nsims, prior_mean, prior_sd) # draw parameters from the prior 
  
  ## Generate data
  if (use_arctanh) {
    phi_samples <- tanh(param_samples)
  } else { # use logit
    phi_samples <- exp(param_samples) / (1 + exp(-param_samples))
  }

  png("./plots/lgss_train_params.png")
  hist(phi_samples, main = "phi")
  dev.off()

  ## Simulate data corresponding to each parameter sample
  data <- lapply(1:nsims, function(i) sim_data(phi = phi_samples[i], iters = iters))

  data_y <- lapply(data, function(x) x$y)
  data_x <- lapply(data, function(x) x$x)
  y_mat <- matrix(unlist(data_y), nrow = nsims, byrow = T)
  x_mat <- matrix(unlist(data_x), nrow = nsims, byrow = T)

  # png("output/lgss_train_data.png")
  # plot(y_mat[90,], type = "l")
  # lines(y_mat[89,], col = "red")
  # dev.off()

  ## Standardise input (time series data) and output (transformed model parameters)
  input_mean <- mean(y_mat)
  input_sd <- sd(y_mat)
  y_mat_std <- (y_mat - input_mean) / input_sd

  ## Split into training, validation and test data (89/10/1)
  inds <- 1:nsims
  train_inds <- sample(inds, 0.89 * length(inds))
  val_inds <- sample(setdiff(inds, train_inds), 0.10 * length(inds))
  test_inds <- setdiff(inds, c(train_inds, val_inds))

  train_input <- y_mat_std[train_inds, ]
  val_input <- y_mat_std[val_inds, ]
  test_input <- y_mat_std[test_inds,]

  train_input <- array_reshape(train_input, c(dim(train_input), 1))
  val_input <- array_reshape(val_input, c(dim(val_input), 1))
  test_input <- array_reshape(test_input, c(dim(test_input), 1))

  train_output <- param_samples[train_inds]
  val_output <- param_samples[val_inds]
  test_output <- param_samples[test_inds]

  simulated_data <- list(train_input = train_input, val_input = val_input, test_input = test_input, 
                        input_mean = input_mean, input_sd = input_sd,
                        train_output = train_output, val_output = val_output, test_output = test_output)

  ## Also need to save the true state for the test data
  train_x_true <- x_mat[train_inds, ]
  val_x_true <- x_mat[val_inds, ]
  test_x_true <- x_mat[test_inds, ]

  true_states <- list(train_x_true = train_x_true, val_x_true = val_x_true, test_x_true = test_x_true)

  if (save_data) {
    qsave(simulated_data, file = paste0("./data/lgss_data_", data_date, ".qs"))
    qsave(true_states, file = paste0("./data/lgss_true_states_", data_date, ".qs"))
  }

} else {
  simulated_data <- qread(paste0("./data/data_", data_date, ".qs"))
  train_input <- simulated_data$train_input
  val_input <- simulated_data$val_input
  test_input <- simulated_data$test_input
  train_output <- simulated_data$train_output
  val_output <- simulated_data$val_output
  test_output <- simulated_data$test_output
}

n_samples <- 3
samples <- sample(1:length(train_inds), n_samples)

png("./plots/lgss_train_data.png", width = 800, height = 600)
par(mfrow = c(n_samples, 1))

for (s in samples) {
  
  train_input_s <- train_input[s,,]
  train_output_s <- train_output[s]
  if (use_arctanh) {
    phi_s <- tanh(train_output_s[1])
  } else {
    phi_s <- exp(train_output_s[1]) / (1 + exp(-train_output_s[1]))
  }
  plot(train_input_s[1:200], type = "l", main = paste("phi = ", round(phi_s, 2)))
}

dev.off()












