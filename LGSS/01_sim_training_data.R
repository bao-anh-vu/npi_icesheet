## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/")

source("./source/sim_data.R")
# source("source/create_model.R")
# source("source/posterior_loss.R")

library(mvtnorm)
# library(keras)
# reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
# library(ggplot2)
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

## Generate data
data_date <- "20250108"

sigma_eps <- 0.5 # measurement error var
sigma_eta <- 0.7 # process error var
phi <- 0.9

# Generate "real" process x_{1:T} and data y_{1:T}
iters <- 500 #T

x_true <- c()
x_true[1] <- 0 #rnorm(1, 0, sqrt(sigma_eta^2 / (1-phi^2)))

for (t in 2:iters) {
  x_true[t] <- phi * x_true[t-1] + rnorm(1, 0, sigma_eta)
}

# Generate observations y_{1:T}
y <- x_true + rnorm(iters, 0, sigma_eps)

## Plot true process and observations
par(mfrow = c(1, 1))
plot(x_true, type = "l")
points(y, col = "cyan")

qsave(y, file = paste0("output/y_", data_date, ".qs"))

## Simulate training data
nsims <- 10000

if (regenerate_data) {
  ## Prior for the parameters: arctanh(phi) ~ N(0, 1), log(sigma_eta) ~ N(0, 1), log(sigma_eps) ~ N(0, 1)
  param_samples <- rmvnorm(nsims, c(0, -1, -1), diag(3)) # draw parameters from the prior 
  
  ## Generate data
  phi_samples <- tanh(param_samples[, 1])
  sigma_eta_samples <- sqrt(exp(param_samples[, 2]))
  sigma_eps_samples <- sqrt(exp(param_samples[, 3]))

  png("output/lgss_train_params.png")
  par(mfrow = c(3, 1))
  hist(phi_samples, main = "phi")
  hist(sigma_eta_samples, main = "sigma_eta")
  hist(sigma_eps_samples, main = "sigma_eps")
  dev.off()

  # 95% credible intervals
  # print(quantile(sigma_eta_samples, c(0.025, 0.975)))

  ## Initial state
  # x_ini_samples <- rnorm(nsims, 0, 1)
  x_ini <- 0 # fix for now

  ## Simulate data corresponding to each parameter sample
  data <- lapply(1:nsims, function(i) sim_data(x_ini = x_ini, 
                                              phi = phi_samples[i], 
                                              sigma_eta = sigma_eta_samples[i], 
                                              sigma_eps = sigma_eps_samples[i], 
                                              iters = iters))
  data_mat <- matrix(unlist(data), nrow = nsims, byrow = T)

  # png("output/lgss_train_data.png")
  # plot(data_mat[90,], type = "l")
  # lines(data_mat[89,], col = "red")
  # dev.off()

  ## Standardise input (time series data) and output (transformed model parameters)
  mean_input <- mean(data_mat)
  sd_input <- sd(data_mat)
  data_mat_std <- (data_mat - mean_input) / sd_input

  ## Output is already standardised since they are draws from N(0, 1)

  ## Split into training, validation and test data (89/10/1)
  inds <- 1:nsims
  train_inds <- sample(inds, 0.89 * length(inds))
  val_inds <- sample(setdiff(inds, train_inds), 0.10 * length(inds))
  test_inds <- setdiff(inds, c(train_inds, val_inds))

  train_input <- data_mat[train_inds, ]
  val_input <- data_mat[val_inds, ]
  test_input <- data_mat[test_inds,]

  train_input <- array_reshape(train_input, c(dim(train_input), 1))
  val_input <- array_reshape(val_input, c(dim(val_input), 1))
  test_input <- array_reshape(test_input, c(dim(test_input), 1))

  train_output <- param_samples[train_inds, ]
  val_output <- param_samples[val_inds, ]
  test_output <- param_samples[test_inds, ]

  # train_output <- cbind(phi_samples[train_inds], sigma_eta_samples[train_inds], sigma_eps_samples[train_inds])
  # val_output <- cbind(phi_samples[val_inds], sigma_eta_samples[val_inds], sigma_eps_samples[val_inds])
  # test_output <- cbind(phi_samples[test_inds], sigma_eta_samples[test_inds], sigma_eps_samples[test_inds])

  simulated_data <- list(train_input = train_input, val_input = val_input, test_input = test_input, 
                          train_output = train_output, val_output = val_output, test_output = test_output)

  if (save_data) {
    qsave(simulated_data, file = paste0("output/lgss_data_", data_date, ".qs"))
  }

} else {
  simulated_data <- qread(paste0("output/data_", data_date, ".qs"))
  train_input <- simulated_data$train_input
  val_input <- simulated_data$val_input
  test_input <- simulated_data$test_input
  train_output <- simulated_data$train_output
  val_output <- simulated_data$val_output
  test_output <- simulated_data$test_output
}

png("output/lgss_train_data.png")
par(mfrow = c(3, 2))

for (i in 15:17) {
  plot(train_input[i, 1:100, ], type = "l")
  plot(phi_samples[i], sigma_eta_samples[i], xlab = "phi", ylab = "sigma_eta", xlim = c(-1, 1), ylim = c(0, 2))
}

dev.off()




