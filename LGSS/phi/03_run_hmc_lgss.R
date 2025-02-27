## Linear Gaussian state space model example
setwd("~/SSA_model/CNN/LGSS/phi/")

library(cmdstanr)
library(qs)

test_on_train <- F
use_arctanh <- T

## Read data
data_date <- "20250108"
simulated_data <- qread(paste0("./data/lgss_data_", data_date, ".qs"))

if (test_on_train) {
  test_input <- simulated_data$train_input#[1:5000, , ]
  test_output <- simulated_data$train_output#[1:5000, ]
} else {
    test_input <- simulated_data$test_input
    test_output <- simulated_data$test_output
}

input_mean <- simulated_data$input_mean
input_sd <- simulated_data$input_sd

# true_sigma_eta <- sqrt(exp(test_output[, 2]))
# true_sigma_eps <- sqrt(exp(test_output[, 3]))

test_sample <- as.numeric(commandArgs(trailingOnly = TRUE))
cat("Test sample: ", test_sample, "\n")

iters <- 10000
burn_in <- 5000
n_chains <- 2
# run_hmc_lgss <- function(data, iters = 10000, burn_in = 5000, n_chains = 1) {
  
y <- test_input[test_sample, , 1] * input_sd + input_mean
stan_file <- "./source/stan_lgss.stan"

lgss_model <- cmdstan_model(
stan_file,
cpp_options = list(stan_threads = TRUE)
)

# log_kappa2_est <- mean(log(y^2)) - (digamma(1/2) + log(2))

lgss_data <- list(Tfin = length(y), y = y, 
                    sigma_eta = 0.7, sigma_eps = 0.5,
                    use_arctanh = ifelse(use_arctanh, 1, 0))

# hfit <- stan(model_code = sv_code, 
#              model_name="sv", data = sv_data, 
#              iter = iters, warmup = burn_in, chains=1)

fit_stan_lgss <- lgss_model$sample(
lgss_data,
chains = n_chains,
parallel_chains = n_chains,
threads = parallel::detectCores(),
refresh = 1000,
iter_warmup = burn_in,
iter_sampling = iters
)

stan_results <- list(draws = fit_stan_lgss$draws(variables = c("phi")),
                    time = fit_stan_lgss$time,
                    summary = fit_stan_lgss$cmdstan_summary,
                    burn_in = burn_in)

## Format test sample to 2 digits
# test_sample <- formatC(test_sample, width = 2, format = "d", flag = "0")

qsave(stan_results, paste0("output/lgss_hmc_results_", 
                            formatC(test_sample, width = 2, format = "d", flag = "0"), 
                            "_", data_date, ".qs"))
#   return(stan_results)
# }

true_vals <- test_output[test_sample]

if (use_arctanh) {
  true_phi <- tanh(true_vals[1])
} else {
  true_phi <- exp(true_vals[1]) / (1 + exp(-true_vals[1]))
}

png(paste0("./plots/lgss_hmc_post_", test_sample, "_", data_date, ".png"), width = 1000, height = 400)

plot(density(stan_results$draws[,,1]), main = "phi")
abline(v = true_phi, col = "black", lty = 2)
dev.off()

png(paste0("./plots/lgss_hmc_trace_", test_sample, "_", data_date, ".png"), width = 1000, height = 400)
plot(stan_results$draws[,,1], type = "l", main = "phi")
abline(h = true_phi, col = "black", lty = 2)
dev.off()
