## Pre-process data

setwd("/home/babv971/SSA_model/CNN/pilot")

rm(list = ls())

save_data <- T 

## Read data
data_date <- "20240320" #"20220329" 
set <- 1
thickness_velocity_arr <- readRDS(file = paste0("./training_data/thickness_velocity_arr_", set, "_", data_date))
friction_arr <- readRDS(file = paste0("./training_data/friction_arr_", set, "_", data_date))
gl_arr <- readRDS(file = paste0("./training_data/gl_arr_", set, "_", data_date))
friction_basis <- readRDS(file = paste0("./training_data/friction_basis_", set, "_", data_date))
basis_coefs <- friction_basis$basis_coefs

## From here on: u = velocity, h = ice thickness

## Standardise output
standard_u_h <- apply(thickness_velocity_arr, 4, function(x) (x - mean(x)) / sd(x), simplify = F)
standard_u_h <- array(unlist(standard_u_h), dim = dim(thickness_velocity_arr))
u_h_mean <- apply(thickness_velocity_arr, 4, mean) # for saving purposes
u_h_sd <- apply(thickness_velocity_arr, 4, sd)

## Now split into train/validation/test sets (89/10/1) 
train_ind <- sample(1:dim(standard_u_h)[1], 0.89*dim(standard_u_h)[1])
val_ind <- sample(setdiff(1:dim(standard_u_h)[1], train_ind), 0.1*dim(standard_u_h)[1])
test_ind <- setdiff(1:dim(standard_u_h)[1], c(train_ind, val_ind))

train_u_h <- standard_u_h[train_ind, , ,]
val_u_h <- standard_u_h[val_ind, , ,]
test_u_h <- standard_u_h[test_ind, , ,]

train_fric <- basis_coefs[train_ind, ]
val_fric <- basis_coefs[val_ind, ]
test_fric <- basis_coefs[test_ind, ]

## Save output
train_data <- list(
    train_input = train_u_h, val_input = val_u_h, test_input = test_u_h,
    input_mean = u_h_mean, input_sd = u_h_sd,
    train_output = train_fric, val_output = val_fric, test_output = test_fric
)

if (save_data) {
    saveRDS(train_data, file = paste0("./training_data/train_data_", data_date))
}

