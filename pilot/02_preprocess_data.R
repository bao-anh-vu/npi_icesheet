## Pre-process data

setwd("/home/babv971/SSA_model/CNN/pilot/")

rm(list = ls())


## Flags
output_var <- "grounding_line" # "friction" # "grounding_line" # "bed_elevation
save_data <- T

library(abind)
# library(tensorflow)
# reticulate::use_condaenv("myenv", required = TRUE)

source("./source/seq_mean_var.R")

## Read data
data_date <- "20240320" #"20220329" 

arg <- commandArgs(trailingOnly = TRUE)
sets <- 1:10
setf <- lapply(sets, function(x) formatC(x, width=2, flag="0"))
# setsf <- paste0("sets", sets[1], "-", sets[length(sets)])#formatC(sets, width=2, flag="0

## Read thickness and velocity data
print("Reading thickness and velocity data...")
t1 <- proc.time()
files <- lapply(setf, function(x) paste0("./training_data/thickness_velocity_arr_", x, "_", data_date, ".rds"))
thickness_velocity_list <- lapply(files, readRDS)
thickness_velocity_arr <- abind(thickness_velocity_list, along = 1)
t2 <- proc.time()

u_mean <- mean(thickness_velocity_list[[1]][,,,1]) #just use the mean from the first set
h_mean <- mean(thickness_velocity_list[[1]][,,,2])
u_sd <- sd(thickness_velocity_list[[1]][,,,1])
h_sd <- sd(thickness_velocity_list[[1]][,,,2])

# u_mean2 <- mean(thickness_velocity_list[[2]][,,,1])
# h_mean2 <- mean(thickness_velocity_list[[2]][,,,2])
# u_sd2 <- sd(thickness_velocity_list[[2]][,,,1])
# h_sd2 <- sd(thickness_velocity_list[[2]][,,,2])

# browser()

rm(thickness_velocity_list)

## Read friction data
print("Reading friction data...")
files <- lapply(setf, function(x) paste0("./training_data/friction_arr_", x, "_", data_date, ".rds"))
fric_list <- lapply(files, readRDS)
fric_arr <- abind(fric_list, along = 1)
rm(fric_list)

## Read grounding line data
print("Reading grounding line data...")
files <- lapply(setf, function(x) paste0("./training_data/gl_arr_", x, "_", data_date, ".rds"))
gl_list <- lapply(files, readRDS)
gl_arr <- abind(gl_list, along = 1)
rm(gl_list)

## Read basis coefficients data
print("Reading basis coefficients data...")
files <- lapply(setf, function(x) paste0("./training_data/friction_basis_", x, "_", data_date, ".rds"))
basis_coefs_list <- lapply(files, function(x) readRDS(x)$basis_coefs)
basis_coefs <- abind(basis_coefs_list, along = 1)
rm(basis_coefs_list)

curr_mem_usage <- sum(sapply(ls(),function(x){object.size(get(x))})) / 1e09


## From here on: u = velocity, h = ice thickness

## Standardise input
print("Standardising input...")

# u_mean <- compute_mean_seq(thickness_velocity_arr[,,,1])
# u_sd <- sqrt(compute_var_seq(thickness_velocity_arr[,,,1]))
std_u <- (thickness_velocity_arr[,,,1] - u_mean) / u_sd

# h_mean <- compute_mean_seq(thickness_velocity_arr[,,,2])
# h_sd <- sqrt(compute_var_seq(thickness_velocity_arr[,,,2]))
std_h <- (thickness_velocity_arr[,,,2] - h_mean) / h_sd

std_input <- abind(std_u, std_h, along = 4)

rm(std_u, std_h)

## Standardise output
print("Standardising output...")

if (output_var == "friction") {
    output <- basis_coefs
} else if (output_var == "grounding_line") {
    output <- drop(gl_arr)
} else if (output_var == "bed_elevation") {
    stop("We can't predict bed elevation yet") 
} else {
    stop("Invalid output_var")
}

# output <- cbind(drop(basis_coefs), drop(gl_arr))

output_mean <- mean(output)
output_sd <- sd(output)
std_output <- (output - output_mean) / output_sd

## Now split into train/validation/test sets (89/10/1) 
print("Splitting data into train/val/test...")
set.seed(2024)
train_ind <- sample(1:dim(std_input)[1], 0.89*dim(std_input)[1])
val_ind <- sample(setdiff(1:dim(std_input)[1], train_ind), 0.1*dim(std_input)[1])
test_ind <- setdiff(1:dim(std_input)[1], c(train_ind, val_ind))

train_input <- std_input[train_ind, , ,]
val_input <- std_input[val_ind, , ,]
test_input <- std_input[test_ind, , ,]
# rm(std_input) # to save memory

# train_fric <- basis_coefs[train_ind, ]
# val_fric <- basis_coefs[val_ind, ]
# test_fric <- basis_coefs[test_ind, ]

train_output <- std_output[train_ind, ]
val_output <- std_output[val_ind, ]
test_output <- std_output[test_ind, ]
# rm(basis_coefs) # save memory

# train_gl <- drop(gl_arr[train_ind, , ,]) # drop dimensions that are just 1, so that we end up with an array of size n_sim x n_years
# val_gl <- drop(gl_arr[val_ind, , ,])
# test_gl <- drop(gl_arr[test_ind, , ,])
#rm(gl_arr) # save memory

## Save output
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])#formatC(sets, width=2, flag="0")

data_dir <- paste0("./training_data/", output_var, "/", setsf) 
if (!dir.exists(data_dir)) {
    dir.create(data_dir)
}

# train_data <- list(
#     train_input = train_input, val_input = val_input, test_input = test_input,
#     input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
#     train_output = train_fric, val_output = val_fric, test_output = test_fric,
#     train_gl = train_gl, val_gl = val_gl, test_gl = test_gl
# )

train_data <- list(
    input = train_input, 
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = train_output,
    output_mean = output_mean, output_sd = output_sd,
    grounding_line = drop(gl_arr[train_ind, , ,])
)

if (save_data) {
    print("Saving output...")
    saveRDS(train_data, file = paste0(data_dir, "/train_data_", data_date, ".rds"))
}
# rm(train_data)

val_data <- list(
    input = val_input, 
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = val_output,
    output_mean = output_mean, output_sd = output_sd,
    grounding_line = drop(gl_arr[val_ind, , ,])
)

if (save_data) {
    print("Saving output...")
    saveRDS(val_data, file = paste0(data_dir, "/val_data_", data_date, ".rds"))
}
# rm(val_data)

test_data <- list(
    input = test_input,
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = test_output,
    output_mean = output_mean, output_sd = output_sd,
    grounding_line = drop(gl_arr[test_ind, , ,])
)

if (save_data) {
    print("Saving output...")
    saveRDS(test_data, file = paste0(data_dir, "/test_data_", data_date, ".rds"))
}

# rm(test_data)
