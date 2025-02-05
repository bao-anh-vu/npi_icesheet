## State inference for the LGSS model conditional on CNN output

source("kalman_filter_1d_lgss.R")

## Read data
data_date <- "20250108"
simulated_data <- qread(paste0("output/lgss_data_", data_date, ".qs"))
train_input <- simulated_data$train_input
val_input <- simulated_data$val_input
test_input <- simulated_data$test_input
train_output <- simulated_data$train_output
val_output <- simulated_data$val_output
test_output <- simulated_data$test_output

