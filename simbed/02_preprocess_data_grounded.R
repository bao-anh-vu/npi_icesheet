## Pre-process data

setwd("~/SSA_model/CNN/simbed/")

rm(list = ls())

## Flags
# sim_beds <- T
# output_var <- "bed" # "friction" # "grounding_line" # "bed_elevation
save_data <- T
standardise_output <- T

library(abind)
# library(tensorflow)
# reticulate::use_condaenv("myenv", required = TRUE)

source("./source/seq_mean_var.R")

## Read data
data_date <- "20240320" # "20220329"

arg <- commandArgs(trailingOnly = TRUE)
sets <- 1#:50 #c(1,3,5) #1:5 #10
setf <- lapply(sets, function(x) formatC(x, width = 2, flag = "0"))
# setsf <- paste0("sets", sets[1], "-", sets[lenhgth(sets)])#formatC(sets, width=2, flag="0

# if (sim_beds) {
    train_data_dir <- "./training_data"
# } else {
    # train_data_dir <- "./training_data"
# }

## Read thickness and velocity data
print("Reading surface data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/surface_obs_arr_", x, "_", data_date, ".rds"))
surface_obs_list <- lapply(files, readRDS)
surface_obs_arr <- abind(surface_obs_list, along = 1)

# first year is the initial condition, so we don't use it as part of the training data
# years <- dim(surface_obs_arr)[3] - 1 
# surface_obs_arr <- surface_obs_arr[, , 2:(years+1), ] 

## Need grounding line positions for each set
## Read grounding line data
print("Reading grounding line data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/gl_arr_", x, "_", data_date, ".rds"))
gl_list <- lapply(files, readRDS)
gl_arr <- abind(gl_list, along = 1)
rm(gl_list)

J <- dim(surface_obs_arr)[2]
gl_points <- floor(gl_arr / 800 * J) # convert from km to grid points
min_gl <- min(gl_points)

## Then subset the surface data up to the grounding line
surface_obs_gr <- surface_obs_arr[, 1:min_gl, , ] # only need the first set of grounding line positions

## might just have to truncate to minimum GL position, otherwise input data has uneven dimensions
## which is not ideal, and might be better to use a mask to ignore the data beyond the grounding line

# compute mean and sd of surface elevation and velocity
if (dim(surface_obs_gr)[1] >= 10000) {
    surf_elev_mean <- mean(surface_obs_gr[1:10000, , , 1]
}
surf_elev_mean <- mean(surface_obs_list[[1]][, , , 1]) # just use the mean from the first set
velocity_mean <- mean(surface_obs_list[[1]][, , , 2])
surf_elev_sd <- sd(surface_obs_list[[1]][, , , 1])
velocity_sd <- sd(surface_obs_list[[1]][, , , 2])

# velocity_mean2 <- mean(surface_obs_list[[2]][,,,1])
# surf_elev_mean2 <- mean(surface_obs_list[[2]][,,,2])
# velocity_sd2 <- sd(surface_obs_list[[2]][,,,1])
# surf_elev_sd2 <- sd(surface_obs_list[[2]][,,,2])

rm(surface_obs_list)

## Read true surface elevation data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_surface_elevs_", x, "_", data_date, ".rds"))
true_surface_list <- lapply(files, readRDS)
true_surface_arr <- abind(true_surface_list, along = 1)
rm(true_surface_list)

## Read true thickness data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_thicknesses_", x, "_", data_date, ".rds"))
true_thickness_list <- lapply(files, readRDS)
true_thickness_arr <- abind(true_thickness_list, along = 1)
rm(true_thickness_list)

## Read true velocity data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_velocities_", x, "_", data_date, ".rds"))
true_velocity_list <- lapply(files, readRDS)
true_velocity_arr <- abind(true_velocity_list, along = 1)
rm(true_velocity_list)

# if (output_var == "friction") {
## Read friction data
print("Reading friction data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/friction_arr_", x, "_", data_date, ".rds"))
fric_list <- lapply(files, readRDS)
fric_arr <- abind(fric_list, along = 1)
rm(fric_list)

## Read basis coefficients data
print("Reading friction basis coefficients data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/friction_basis_", x, "_", data_date, ".rds"))
fric_basis_coefs_list <- lapply(files, function(x) readRDS(x)$basis_coefs)
fric_basis_coefs <- abind(fric_basis_coefs_list, along = 1)
fric_basis_mat <- readRDS(files[[1]])$basis_mat
rm(fric_basis_coefs_list)

# } else if (output_var == "bed") {
# Read true bed data
print("Reading bed data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/bed_arr_", x, "_", data_date, ".rds"))
bed_list <- lapply(files, readRDS)
bed_arr <- abind(bed_list, along = 1)
rm(bed_list)

## Read basis coefficients data
print("Reading bed basis coefficients data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/bed_basis_", x, "_", data_date, ".rds"))
bed_basis_coefs_list <- lapply(files, function(x) readRDS(x)$basis_coefs)
bed_basis_coefs <- abind(bed_basis_coefs_list, along = 1)
bed_basis_mat <- readRDS(files[[1]])$basis_mat
rm(bed_basis_coefs_list)

# }

curr_mem_usage <- sum(sapply(ls(), function(x) {
    object.size(get(x))
})) / 1e09

## From here on: u = velocity, h = ice thickness

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
# print("Standardising output...")

# if (output_var == "friction") {
#     output <- basis_coefs
# } else if (output_var == "grounding_line") {
#     output <- drop(gl_arr)
# } else if (output_var == "bed") {
#     output <- basis_coefs
# } else {
#     stop("Invalid output_var")
# }

# output_mean <- mean(output)
# output_sd <- sd(output)
# std_output <- (output - output_mean) / output_sd

# fric_mean <- mean(fric_basis_coefs)
# fric_sd <- sd(fric_basis_coefs)
# std_fric <- (fric_basis_coefs - fric_mean) / fric_sd

# fric_mean <- mean(fric_basis_coefs)
# fric_sd <- sd(fric_basis_coefs)
# std_fric <- (fric_basis_coefs - fric_mean) / fric_sd

## Now split into train/validation/test sets (89/10/1)
print("Splitting data into train/val/test...")
set.seed(2024)
train_ind <- sample(1:dim(std_input)[1], 0.89 * dim(std_input)[1])
val_ind <- sample(setdiff(1:dim(std_input)[1], train_ind), 0.1 * dim(std_input)[1])
test_ind <- setdiff(1:dim(std_input)[1], c(train_ind, val_ind))

train_input <- std_input[train_ind, , , ]
val_input <- std_input[val_ind, , , ]
test_input <- std_input[test_ind, , , ]
# rm(std_input) # to save memory

# train_fric <- fric_basis_coefs[train_ind, ]
# val_fric <- fric_basis_coefs[val_ind, ]
# test_fric <- fric_basis_cofefs[test_ind, ]

# train_bed <- bed_basis_coefs[train_ind, ]
# val_bed <- bed_basis_coefs[val_ind, ]
# test_bed <- bed_basis_coefs[test_ind, ]


# train_output <- std_output[train_ind, ]
# val_output <- std_output[val_ind, ]
# test_output <- std_output[test_ind, ]
# rm(basis_coefs) # save memory

# train_gl <- drop(gl_arr[train_ind, , ,]) # drop dimensions that are just 1, so that we end up with an array of size n_sim x n_years
# val_gl <- drop(gl_arr[val_ind, , ,])
# test_gl <- drop(gl_arr[test_ind, , ,])
# rm(gl_arr) # save memory


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

## Save output
setsf <- paste0("sets", sets[1], "-", sets[length(sets)]) # formatC(sets, width=2, flag="0")

data_dir <- paste0(train_data_dir, "/", setsf)

if (!dir.exists(data_dir)) {
    dir.create(data_dir)
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
    true_surface_elev = true_surface_arr[train_ind, , ,],
    true_thickness_train =  true_thickness_arr[train_ind, , ,],
    true_velocity_train = true_velocity_arr[train_ind, , ,]
)

if (save_data) {
    print("Saving training data...")
    saveRDS(train_data, file = paste0(data_dir, "/train_data_", data_date, ".rds"))
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
    true_surface_elev = true_surface_arr[val_ind, , ,],
    true_thickness_val =  true_thickness_arr[val_ind, , ,],
    true_velocity_val = true_velocity_arr[val_ind, , ,]
)

if (save_data) {
    print("Saving validation data...")
    saveRDS(val_data, file = paste0(data_dir, "/val_data_", data_date, ".rds"))
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
    true_surface_elev = true_surface_arr[test_ind, , ,],
    true_thickness_test =  true_thickness_arr[test_ind, , ,],
    true_velocity_test = true_velocity_arr[test_ind, , ,]
)

if (save_data) {
    print("Saving test data...")
    saveRDS(test_data, file = paste0(data_dir, "/test_data_", data_date, ".rds"))
}

# rm(test_data)
