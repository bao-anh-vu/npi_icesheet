## Pre-process data

setwd("/home/babv971/SSA_model/CNN/pilot/")

rm(list = ls())

save_data <- T

library(abind)
# library(tensorflow)
# reticulate::use_condaenv("myenv", required = TRUE)

source("./source/seq_mean_var.R")


## Read data
data_date <- "20240320" #"20220329" 

# arg <- commandArgs(trailingOnly = TRUE)
sets <- 1
setf <- lapply(sets, function(x) formatC(x, width=2, flag="0"))
# setsf <- paste0("sets", sets[1], "-", sets[length(sets)])#formatC(sets, width=2, flag="0

## Read thickness and velocity data
print("Reading thickness and velocity data...")
t1 <- proc.time()
files <- lapply(setf, function(x) paste0("./training_data/thickness_velocity_arr_", x, "_", data_date, ".rds"))
thickness_velocity_list <- lapply(files, readRDS)
thickness_velocity_arr <- abind(thickness_velocity_list, along = 1)
t2 <- proc.time()

u_mean <- mean(thickness_velocity_list[[1]][,,,1])
h_mean <- mean(thickness_velocity_list[[1]][,,,2])
u_sd <- sd(thickness_velocity_list[[1]][,,,1])
h_sd <- sd(thickness_velocity_list[[1]][,,,2])
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

# setf <- formatC(sets[1], width=2, flag="0")
# thickness_velocity_arr <- readRDS(file = paste0("./training_data/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
# friction_arr <- readRDS(file = paste0("./training_data/friction_arr_", setf, "_", data_date, ".rds"))
# gl_arr <- readRDS(file = paste0("./training_data/gl_arr_", setf, "_", data_date, ".rds"))
# fric_basis <- readRDS(file = paste0("./training_data/friction_basis_", setf, "_", data_date, ".rds"))
# basis_coefs <- fric_basis$basis_coefs
# for (set in sets[2:length(sets)]) {
#     cat("Set", set, "\n")
#     setf <- formatC(set, width=2, flag="0")
#     u_h_set <- readRDS(file = paste0("./training_data/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
#     fric_set <- readRDS(file = paste0("./training_data/friction_arr_", setf, "_", data_date, ".rds"))
#     gl_set <- readRDS(file = paste0("./training_data/gl_arr_", setf, "_", data_date, ".rds"))
#     fric_basis_set <- readRDS(file = paste0("./training_data/friction_basis_", setf, "_", data_date, ".rds"))
#     basis_coefs_set <- fric_basis_set$basis_coefs

#     thickness_velocity_arr <- abind(thickness_velocity_arr, u_h_set, along = 1)
#     friction_arr <- rbind(friction_arr, fric_set, along = 1)
#     gl_arr <- abind(gl_arr, gl_set, along = 1)
#     basis_coefs <- rbind(basis_coefs, basis_coefs_set)
# }

curr_mem_usage <- sum(sapply(ls(),function(x){object.size(get(x))})) / 1e09


# saveRDS(thickness_velocity_arr, file = paste0("./training_data/thickness_velocity_arr_", setsf, "_", data_date, ".rds"))
# thickness_velocity_arr <- readRDS(file = paste0("./training_data/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
# friction_arr <- readRDS(file = paste0("./training_data/friction_arr_", setf, "_", data_date, ".rds"))
# gl_arr <- readRDS(file = paste0("./training_data/gl_arr_", setf, "_", data_date, ".rds"))
# friction_basis <- readRDS(file = paste0("./training_data/friction_basis_", setf, "_", data_date, ".rds"))
# basis_coefs <- friction_basis$basis_coefs

## From here on: u = velocity, h = ice thickness

## Standardise output
print("Standardising data...")

# u_mean <- compute_mean_seq(thickness_velocity_arr[,,,1])
# u_sd <- sqrt(compute_var_seq(thickness_velocity_arr[,,,1]))
standard_u <- (thickness_velocity_arr[,,,1] - u_mean) / u_sd

# h_mean <- compute_mean_seq(thickness_velocity_arr[,,,2])
# h_sd <- sqrt(compute_var_seq(thickness_velocity_arr[,,,2]))
standard_h <- (thickness_velocity_arr[,,,2] - h_mean) / h_sd

standard_u_h <- abind(standard_u, standard_h, along = 4)

rm(standard_u, standard_h)
# mean_10k <- mean(thickness_velocity_arr[1:10000,,,1])
# sd_10k <- sd(thickness_velocity_arr[1:10000,,,1])

# browser()

# standard_u_h <- apply(thickness_velocity_arr, 1, scale, simplify = F)


# standard_u_h <- apply(thickness_velocity_arr, 4, 
#                         function(x) (x - mean(x)) / sd(x), simplify = F)
# standard_u_h <- array(unlist(standard_u_h), dim = dim(thickness_velocity_arr))

## Now split into train/validation/test sets (89/10/1) 
print("Splitting data into train/val/test...")
set.seed(2024)
train_ind <- sample(1:dim(standard_u_h)[1], 0.89*dim(standard_u_h)[1])
val_ind <- sample(setdiff(1:dim(standard_u_h)[1], train_ind), 0.1*dim(standard_u_h)[1])
test_ind <- setdiff(1:dim(standard_u_h)[1], c(train_ind, val_ind))

train_u_h <- standard_u_h[train_ind, , ,]
val_u_h <- standard_u_h[val_ind, , ,]
test_u_h <- standard_u_h[test_ind, , ,]
rm(standard_u_h) # to save memory

train_fric <- basis_coefs[train_ind, ]
val_fric <- basis_coefs[val_ind, ]
test_fric <- basis_coefs[test_ind, ]
rm(basis_coefs) # save memory

train_gl <- as.vector(gl_arr[train_ind, , ,])
val_gl <- as.vector(gl_arr[val_ind, , ,])
test_gl <- as.vector(gl_arr[test_ind, , ,])
rm(gl_arr) # save memory

## Save output
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])#formatC(sets, width=2, flag="0")

train_data_path <- paste0("training_data/", setsf)

if (!file.exists(train_data_path)) {
    dir.create(train_data_path)
}

# train_data <- list(
#     train_input = train_u_h, val_input = val_u_h, test_input = test_u_h,
#     input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
#     train_output = train_fric, val_output = val_fric, test_output = test_fric,
#     train_gl = train_gl, val_gl = val_gl, test_gl = test_gl
# )

train_data <- list(
    input = train_u_h, 
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = train_fric, 
    grounding_line = train_gl
)

if (save_data) {
    print("Saving output...")
    saveRDS(train_data, file = paste0(train_data_path, "/train_data_", data_date, ".rds"))
}
rm(train_data)

val_data <- list(
    input = val_u_h, 
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = val_fric, 
    grounding_line = val_gl
)

if (save_data) {
    print("Saving output...")
    saveRDS(val_data, file = paste0(train_data_path, "/val_data_", data_date, ".rds"))
}
rm(val_data)

test_data <- list(
    input = test_u_h,
    input_mean = c(u_mean, h_mean), input_sd = c(u_sd, h_sd),
    output = test_fric,
    grounding_line = test_gl
)

if (save_data) {
    print("Saving output...")
    saveRDS(test_data, file = paste0(train_data_path, "/test_data_", data_date, ".rds"))
}

rm(test_data)
