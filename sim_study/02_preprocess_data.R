## Pre-process data

setwd("~/SSA_model/CNN/simbed/")

rm(list = ls())

# library(keras)
# reticulate::use_condaenv("myenv", required = TRUE)
# library(tensorflow)
library(ggplot2)
library(dplyr)
library(abind)
library(parallel)
library(qs)

## Flags
# sim_beds <- T
# output_var <- "bed" # "friction" # "grounding_line" # "bed_elevation
save_data <- T
standardise_output <- T
use_missing_pattern <- T

# library(tensorflow)
# reticulate::use_condaenv("myenv", required = TRUE)

source("./source/seq_mean_var.R")

# # List physical devices
# gpus <- tf$config$experimental$list_physical_devices('GPU')

# if (length(gpus) > 0) {
#   tryCatch({
#     # Restrict TensorFlow to only allocate 4GB of memory on the first GPU
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

## Read data
data_date <- "20240320" # "20220329"

arg <- commandArgs(trailingOnly = TRUE)
sets <- 1:10 #1:50 
setf <- lapply(sets, function(x) formatC(x, width = 2, flag = "0"))
# setsf <- paste0("sets", sets[1], "-", sets[lenhgth(sets)])#formatC(sets, width=2, flag="0

train_data_dir <- "./training_data"

## Read thickness and velocity data
print("Reading surface data...")
files <- lapply(setf, function(x) paste0(train_data_dir, "/surface_obs_arr_", x, "_", data_date, ".qs"))
surface_obs_list <- mclapply(files, qread, mc.cores = 10L)
surface_obs_arr <- abind(surface_obs_list, along = 1)

if (use_missing_pattern) {
    print("Reading missing patterns...")
    surf_elev_missing_pattern <- readRDS("./training_data/surf_elev_missing_pattern.rds")
    vel_missing_pattern <- readRDS("./training_data/vel_missing_pattern.rds")
    missing_patterns <- abind(list(surf_elev_missing_pattern, vel_missing_pattern), along = 3)

    surface_obs_list <- lapply(1:dim(surface_obs_arr)[1], function(i) { surface_obs_arr[i,,,]})

    surface_obs_list_missing <- lapply(surface_obs_list, function(arr) {
    se <- arr[,,1] * surf_elev_missing_pattern
    vel <- arr[,,2] * vel_missing_pattern
    abind(se, vel, along = 3)
    })

    system.time({
        surface_obs_arr <- abind(surface_obs_list_missing, along = 0)
    })

    ## Plot missing pattern
    se_mp_c <- c(surf_elev_missing_pattern)
    se_data_pattern <- c(surface_obs_list_missing[[1]][,,1])
    J <- dim(surface_obs_arr)[2]
    se_mp_df <- data.frame(gridpt = rep(1:J, ncol(surf_elev_missing_pattern)), 
                            nonmissing = se_mp_c,
                            val = se_data_pattern, # apply missing pattern to dataset
                            year = rep(0:(ncol(surf_elev_missing_pattern)-1), each = J))

    se_mp_df <- se_mp_df %>% mutate(val_mp = ifelse(nonmissing == 1, val, NA))

    ## Plot missing surface elevation pattern over space-time
    se_st_plot <- ggplot(se_mp_df) + 
        geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) +
        # geom_tile(aes(x = gridpt, y = year, fill = val_mp)) +
        scale_y_reverse(
            breaks = seq(min(se_mp_df$year)+1, max(se_mp_df$year), 5),
            labels = seq(min(se_mp_df$year)+1, max(se_mp_df$year), 5)
        ) +
        scale_fill_manual(values=c("grey", "turquoise")) +
        theme_bw() + 
        labs(x = "Grid Point", y = "Year", fill = "Non-missing") +
        # labs(x = "Grid Point", y = "Year", fill = "Surface elev.") +
        # scale_fill_distiller(palette = "Blues", direction = 1) +
        # ggtitle("Thwaites Glacier Surface Elevation Over Time") + 
        theme(text = element_text(hjust = 0.5, size = 24))

    png("./plots/surf_elev_missing_pattern.png", width = 1200, height = 800, res = 150)
    print(se_st_plot)
    dev.off()

    vel_mp_c <- c(vel_missing_pattern)
    vel_data_pattern <- c(surface_obs_list_missing[[1]][,,2])
    vel_mp_df <- data.frame(gridpt = rep(1:J, ncol(vel_missing_pattern)), 
                            nonmissing = vel_mp_c, 
                            val = vel_data_pattern,
                            year = rep(0:(ncol(vel_missing_pattern)-1), each = J))

    vel_mp_df <- vel_mp_df %>% mutate(val_mp = ifelse(nonmissing == 1, val, NA))

    ## Plot missing velocity pattern over space-time
    vel_st_plot <- ggplot(vel_mp_df) + 
        geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) + 
        scale_fill_manual(values=c("grey", "turquoise")) +
        scale_y_reverse(
            breaks = seq(min(se_mp_df$year)+1, max(se_mp_df$year), 5),
            labels = seq(min(se_mp_df$year)+1, max(se_mp_df$year), 5)
        ) +
        # geom_tile(aes(x = gridpt, y = year, fill = val_mp)) +
        theme_bw() + 
        labs(x = "Grid Point", y = "Year", fill = "Non-missing") + 
        # labs(x = "Grid Point", y = "Year", fill = "Surface velocity") + 
        # scale_fill_distiller(palette = "Reds", direction = 1) +
        # ggtitle("Thwaites Glacier Velocity Over Time") + 
        theme(text = element_text(hjust = 0.5, size = 24))

    png("./plots/vel_missing_pattern.png", width = 1200, height = 800, res = 150)
    print(vel_st_plot)
    dev.off()

    browser()
}

# compute mean and sd of surface elevation and velocity
if (dim(surface_obs_arr)[1] <= 10000) {
    surf_elev_mean <- mean(surface_obs_arr[, , , 1]) # just use the mean from the first set
    velocity_mean <- mean(surface_obs_arr[, , , 2])
    surf_elev_sd <- sd(surface_obs_arr[, , , 1])
    velocity_sd <- sd(surface_obs_arr[, , , 2])
} else {
    surf_elev_mean <- mean(surface_obs_arr[1:10000, , , 1]) # just use the mean from the first set
    velocity_mean <- mean(surface_obs_arr[1:10000, , , 2])
    surf_elev_sd <- sd(surface_obs_arr[1:10000, , , 1])
    velocity_sd <- sd(surface_obs_arr[1:10000, , , 2])
}

rm(surface_obs_list)

print("Reading ground truth...")
## Read true surface elevation data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_surface_elevs_", x, "_", data_date, ".qs"))
t1 <- system.time({true_surface_list <- mclapply(files, qread, mc.cores = 10L)})

true_surface_arr <- abind(true_surface_list, along = 1)
rm(true_surface_list)

## Read true thickness data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_thicknesses_", x, "_", data_date, ".qs"))
true_thickness_list <- mclapply(files, qread, mc.cores = 10L)
true_thickness_arr <- abind(true_thickness_list, along = 1)
rm(true_thickness_list)

## Read true velocity data
files <- lapply(setf, function(x) paste0(train_data_dir, "/true_velocities_", x, "_", data_date, ".qs"))
true_velocity_list <- mclapply(files, qread, mc.cores = 10L)
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

# }

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

if (use_missing_pattern) {
    data_dir <- paste0(data_dir, "/missing")
} else {
    data_dir <- paste0(data_dir, "/nonmissing")
}

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
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
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
    qsave(train_data, file = paste0(data_dir, "/train_data_", data_date, ".qs"))
}
# rm(train_data)

val_data <- list(
    input = val_input,
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
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
    qsave(val_data, file = paste0(data_dir, "/val_data_", data_date, ".qs"))
}
# rm(val_data)

test_data <- list(
    input = test_input,
    input_mean = c(surf_elev_mean, velocity_mean), input_sd = c(surf_elev_sd, velocity_sd),
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
    qsave(test_data, file = paste0(data_dir, "/test_data_", data_date, ".qs"))
}

# rm(test_data)
