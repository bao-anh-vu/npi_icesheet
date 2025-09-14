## Pre-process data

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(abind)
library(parallel)
library(qs)

## Flags
# sim_beds <- T
# output_var <- "bed" # "friction" # "grounding_line" # "bed_elevation
save_data <- T
standardise_output <- T
use_missing_pattern <- F

# library(tensorflow)
# reticulate::use_condaenv("myenv", required = TRUE)

source("./source/seq_mean_var.R")

# List physical devices
gpus <- tf$config$experimental$list_physical_devices('GPU')

if (length(gpus) > 0) {
  tryCatch({
    # Restrict TensorFlow to only allocate 4GB of memory on the first GPU
    tf$config$experimental$set_virtual_device_configuration(
      gpus[[1]],
      list(tf$config$experimental$VirtualDeviceConfiguration(memory_limit=4096*10))
    )
    
    logical_gpus <- tf$config$experimental$list_logical_devices('GPU')
    
    print(paste0(length(gpus), " Physical GPUs,", length(logical_gpus), " Logical GPUs"))
  }, error = function(e) {
    # Virtual devices must be set before GPUs have been initialized
    print(e)
  })
}

## Read data
data_date <- "20241111" #"20241103" # "20220329"

arg <- commandArgs(trailingOnly = TRUE)
sets <- 1:20 #c(1,3,5) #1:5 #10
setf <- lapply(sets, function(x) formatC(x, width = 2, flag = "0"))
# setsf <- paste0("sets", sets[1], "-", sets[lenhgth(sets)])#formatC(sets, width=2, flag="0

train_data_dir <- "./data/training_data"

## Read surface elevation and velocity hdata
print("Reading surface data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/surface_obs_arr_", x, "_", data_date, ".qs"))
surface_obs_list <- lapply(files, qread)
surface_obs_arr <- abind(surface_obs_list, along = 1)

if (use_missing_pattern) {
    print("Reading missing patterns...")
    surf_elev_missing_pattern <- qread("./data/surface_elev/missing_pattern.qs")
    vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")
    missing_patterns <- abind(list(surf_elev_missing_pattern, vel_missing_pattern), along = 3)

    surface_obs_list <- lapply(1:dim(surface_obs_arr)[1], function(i) { surface_obs_arr[i,,,]})

    ## Multiply the missing patterns by the surface_obs_list
    surface_obs_list_missing <- lapply(surface_obs_list, function(arr) {
    se <- arr[,,1] * surf_elev_missing_pattern
    vel <- arr[,,2] * vel_missing_pattern
    abind(se, vel, along = 3)
    })

    system.time({
        surface_obs_arr <- abind(surface_obs_list_missing, along = 0)
    })
    
}


test <- sapply(1:dim(surface_obs_arr)[1], function(i) ifelse(max(surface_obs_arr[i,,,]) > 10000 | min(surface_obs_arr[i,,,]) < -100, 1, 0))
bad_sims <- which(test == 1)

if (length(bad_sims) > 0) { surface_obs_arr <- surface_obs_arr[-bad_sims, , ,] }

## Calculate the mean of surface_obs_arr over the first 3 dimensions
if (dim(surface_obs_arr)[1] <= 10000) { # if fewer than 10000 simulations, use all of them
    n_for_mean <- 1:dim(surface_obs_arr)[1]
} else { # otherwise, sample 10000 simulations
    n_for_mean <- sample(1:dim(surface_obs_arr)[1], 10000)
}

surf_elev_mean <- mean(surface_obs_arr[n_for_mean, , , 1]) # just use the mean from the first set
velocity_mean <- mean(surface_obs_arr[n_for_mean, , , 2])
surf_elev_sd <- sd(surface_obs_arr[n_for_mean, , , 1])
velocity_sd <- sd(surface_obs_arr[n_for_mean, , , 2])

rm(surface_obs_list)

## Read true surface elevation data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_surface_elevs_", x, "_", data_date, ".qs"))
true_surface_list <- lapply(files, qread)
true_surface_arr <- abind(true_surface_list, along = 1)
rm(true_surface_list)

## Read true thickness data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_thicknesses_", x, "_", data_date, ".qs"))
true_thickness_list <- lapply(files, qread)
true_thickness_arr <- abind(true_thickness_list, along = 1)
rm(true_thickness_list)

## Read true velocity data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_velocities_", x, "_", data_date, ".qs"))
true_velocity_list <- lapply(files, qread)
true_velocity_arr <- abind(true_velocity_list, along = 1)
rm(true_velocity_list)

# if (output_var == "friction") {
## Read friction data
print("Reading friction data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/friction_arr_", x, "_", data_date, ".qs"))
fric_list <- lapply(files, qread)
fric_arr <- abind(fric_list, along = 1)
rm(fric_list)

## Read basis coefficients data
print("Reading friction basis coefficients data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/friction_basis_", x, "_", data_date, ".qs"))
fric_basis_coefs_list <- lapply(files, function(x) qread(x)$basis_coefs)
fric_basis_coefs <- abind(fric_basis_coefs_list, along = 1)
fric_basis_mat <- qread(files[[1]])$basis_mat
rm(fric_basis_coefs_list)

# } else if (output_var == "bed") {
# Read true bed data
print("Reading bed data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/bed_arr_", x, "_", data_date, ".qs"))
bed_list <- lapply(files, qread)
bed_arr <- abind(bed_list, along = 1)
rm(bed_list)

## Read basis coefficients data
print("Reading bed basis coefficients data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/bed_basis_", x, "_", data_date, ".qs"))
bed_basis_coefs_list <- lapply(files, function(x) qread(x)$basis_coefs)
bed_basis_coefs <- abind(bed_basis_coefs_list, along = 1)
bed_basis_mat <- qread(files[[1]])$basis_mat
rm(bed_basis_coefs_list)

## Read grounding line data
print("Reading grounding line data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/gl_arr_", x, "_", data_date, ".qs"))
gl_list <- lapply(files, qread)
gl_arr <- abind(gl_list, along = 1)
rm(gl_list)

curr_mem_usage <- sum(sapply(ls(), function(x) {
    object.size(get(x))
})) / 1e09

## From here on: u = velocity, h = ice thickness

## Filter out bad simulations
if (length(bad_sims) > 0) {
    print(paste0("Removing ", length(bad_sims), " bad simulations..."))
    true_surface_arr <- true_surface_arr[-bad_sims, , ,]
    true_thickness_arr <- true_thickness_arr[-bad_sims, , ,]
    true_velocity_arr <- true_velocity_arr[-bad_sims, , ,]
    fric_arr <- fric_arr[-bad_sims,]
    fric_basis_coefs <- fric_basis_coefs[-bad_sims, ]
    bed_arr <- bed_arr[-bad_sims,]
    bed_basis_coefs <- bed_basis_coefs[-bad_sims, ]
    gl_arr <- gl_arr[-bad_sims, ]
}

## Standardise input
print("Standardising input...")

# velocity_mean <- compute_mean_seq(surface_obs_arr[,,,1])
# velocity_sd <- sqrt(compute_var_seq(surface_obs_arr[,,,1]))
std_z <- (surface_obs_arr[, , , 1] - surf_elev_mean) / surf_elev_sd

# surf_elev_mean <- compute_mean_seq(surface_obs_arr[,,,2])
# surf_elev_sd <- sqrt(compute_var_seq(surface_obs_arr[,,,2]))
std_u <- (surface_obs_arr[, , , 2] - velocity_mean) / velocity_sd 

std_input <- abind(std_z, std_u, along = 4)


rm(std_z, std_u)

## Standardise output
if (standardise_output) {
    # Need to standardise the bed and friction coefs separately here, then combine into one big array
    mean_fric_coefs <- mean(fric_basis_coefs)
    sd_fric_coefs <- sd(fric_basis_coefs)

    mean_bed_coefs <- mean(bed_basis_coefs)
    sd_bed_coefs <- sd(bed_basis_coefs)

    mean_gl <- mean(gl_arr)
    sd_gl <- sd(gl_arr)

    fric_basis_coefs <- (fric_basis_coefs - mean_fric_coefs) / sd_fric_coefs
    bed_basis_coefs <- (bed_basis_coefs - mean_bed_coefs) / sd_bed_coefs
    gl_arr_std <- (gl_arr - mean_gl) / sd_gl

} else {
    mean_fric_coefs <- 0
    sd_fric_coefs <- 1
    nean_bed_coefs <- 0
    sd_bed_coefs <- 1
    mean_gl <- 0
    sd_gl <- 1
}

## Now split into train/validation/test sets (89/10/1)
print("Splitting data into train/val/test...")
set.seed(2024)
remaining_sims <- dim(std_input)[1] - length(bad_sims)
train_ind <- sample(1:remaining_sims, 0.89 * remaining_sims)
val_ind <- sample(setdiff(1:remaining_sims, train_ind), 0.1 * remaining_sims)
test_ind <- setdiff(1:remaining_sims, c(train_ind, val_ind))

train_input <- std_input[train_ind, , , ]
val_input <- std_input[val_ind, , , ]
test_input <- std_input[test_ind, , , ]

## Save output
setsf <- paste0("sets", sets[1], "-", sets[length(sets)]) # formatC(sets, width=2, flag="0")

# if (use_missing_pattern) {
    output_dir <- paste0(train_data_dir, "/", setsf)
# } else {
#     data_dir <- paste0(train_data_dir, "/", setsf)
# }

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

if (use_missing_pattern) {
    output_dir <- paste0(output_dir, "/missing")
} else {
    output_dir <- paste0(output_dir, "/nonmissing")
}

if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

## True friction/bed for each set
true_fric_train <- fric_arr[train_ind, ]
true_fric_val <- fric_arr[val_ind, ]
true_fric_test <- fric_arr[test_ind, ]

true_bed_train <- bed_arr[train_ind, ]
true_bed_val <- bed_arr[val_ind, ]
true_bed_test <- bed_arr[test_ind, ]

true_gl_train <- drop(gl_arr[train_ind, ])
true_gl_val <- drop(gl_arr[val_ind, ])
true_gl_test <- drop(gl_arr[test_ind, ])

train_data <- list(
    input = train_input,
    # input_mean = c(velocity_mean, surf_elev_mean), input_sd = c(velocity_sd, surf_elev_sd),
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
    # output = train_output,
    # output_mean = output_mean, output_sd = output_sd,
    # truth = train_truth,
    bed_coefs = bed_basis_coefs[train_ind, ],
    mean_bed_coefs = mean_bed_coefs,
    sd_bed_coefs = sd_bed_coefs,
    fric_coefs = fric_basis_coefs[train_ind, ],
    mean_fric_coefs = mean_fric_coefs,
    sd_fric_coefs = sd_fric_coefs,
    grounding_line = drop(gl_arr_std[train_ind, ]),
    mean_gl = mean_gl,
    sd_gl = sd_gl,
    true_bed = true_bed_train,
    true_fric = true_fric_train,
    true_gl = true_gl_train,
    true_surface_elev = true_surface_arr[train_ind, , ],
    true_thickness_train =  true_thickness_arr[train_ind, , ],
    true_velocity_train = true_velocity_arr[train_ind, , ]
)

if (save_data) {
    print("Saving training data...")
    qsave(train_data, file = paste0(output_dir, "/train_data_", data_date, ".qs"))
}
# rm(train_data)

val_data <- list(
    input = val_input,
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
    # output = val_output,
    # output_mean = output_mean, output_sd = output_sd,
    # truth = val_truth,
    bed_coefs = bed_basis_coefs[val_ind, ],
    mean_bed_coefs = mean_bed_coefs,
    sd_bed_coefs = sd_bed_coefs,
    fric_coefs = fric_basis_coefs[val_ind, ],
    mean_fric_coefs = mean_fric_coefs,
    sd_fric_coefs = sd_fric_coefs,
    grounding_line = drop(gl_arr_std[val_ind, ]),
    mean_gl = mean_gl,
    sd_gl = sd_gl,
    true_bed = true_bed_val,
    true_fric = true_fric_val,
    true_gl = true_gl_val,
    true_surface_elev = true_surface_arr[val_ind, , ],
    true_thickness_val =  true_thickness_arr[val_ind, , ],
    true_velocity_val = true_velocity_arr[val_ind, , ]
)

if (save_data) {
    print("Saving validation data...")
    qsave(val_data, file = paste0(output_dir, "/val_data_", data_date, ".qs"))
}
# rm(val_data)

test_data <- list(
    input = test_input,
    # input_mean = c(velocity_mean, surf_elev_mean), input_sd = c(velocity_sd, surf_elev_sd),
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
    # output = test_output,
    # output_mean = output_mean, output_sd = output_sd,
    # truth = test_truth,
    bed_coefs = bed_basis_coefs[test_ind, ],
    mean_bed_coefs = mean_bed_coefs,
    sd_bed_coefs = sd_bed_coefs,
    fric_coefs = fric_basis_coefs[test_ind, ],
    mean_fric_coefs = mean_fric_coefs,
    sd_fric_coefs = sd_fric_coefs,
    
    bed_basis_mat = bed_basis_mat,
    fric_basis_mat = fric_basis_mat,
    grounding_line = drop(gl_arr_std[test_ind, ]),
    mean_gl = mean_gl,
    sd_gl = sd_gl,
    true_bed = true_bed_test,
    true_fric = true_fric_test,
    true_gl = true_gl_test,
    true_surface_elev = true_surface_arr[test_ind, , ],
    true_thickness_test =  true_thickness_arr[test_ind, , ],
    true_velocity_test = true_velocity_arr[test_ind, , ]
)

if (save_data) {
    print("Saving test data...")
    qsave(test_data, file = paste0(output_dir, "/test_data_", data_date, ".qs"))
}

# rm(test_data)
