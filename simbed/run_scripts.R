# Rscript 01_simulate_data.R
# setwd("~/SSA_model/CNN/simbed/")

# print("Running script 01_simulate_data...")
# source("01_simulate_data.R")

print("Running script 02_preprocess_data...")
source("02_preprocess_data.R")

print("Running script 03_run_cnn...")
source("03_run_cnn.R")

print("Running script 04_run_cnn_quantiles.t..")
source("04_run_cnn_quantiles.R") 

# print("Running script 05...")
# source("05_post_process_quantiles.R") 

## Note to self: rerunning the case where log_transform = F