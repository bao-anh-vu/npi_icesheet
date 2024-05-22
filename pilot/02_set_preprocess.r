setwd("/home/babv971/SSA_model/CNN/pilot/")

rm(list = ls())

library(abind)

save_data <- T

## Read data
data_date <- "20240320" # "20220329"

sets <- 1:2
setf <- lapply(sets, function(x) formatC(x, width = 2, flag = "0"))

first_set <- sets[1]

## Read thickness and velocity data
print("Reading thickness and velocity data...")

thickness_velocity_arr <- readRDS(paste0(
    "./training_data/thickness_velocity_arr_",
    formatC(first_set, width = 2, flag = "0"), "_", data_date, ".rds"
))

u_h_mean <- apply(thickness_velocity_arr, 4, mean)
u_h_sd <- apply(thickness_velocity_arr, 4, sd)

print("Standardising first set...")
# u_mean <- compute_mean_seq(thickness_velocity_arr[,,,1])
# u_sd <- sqrt(compute_var_seq(thickness_velocity_arr[,,,1]))
standard_u_h <- (thickness_velocity_arr - u_h_mean) / u_h_sd

## Read friction data
print("Reading friction data...")
fric_arr <- readRDS(paste0("./training_data/friction_arr_", formatC(first_set, width = 2, flag = "0"), "_", data_date, ".rds"))

## Read grounding line data
print("Reading grounding line data...")
gl_arr <- readRDS(paste0("./training_data/gl_arr_", formatC(first_set, width = 2, flag = "0"), "_", data_date, ".rds"))

## Read basis coefficients data
print("Reading basis coefficients data...")
fric_basis <- readRDS(paste0("./training_data/friction_basis_", formatC(first_set, width = 2, flag = "0"), "_", data_date, ".rds"))
basis_coefs <- fric_basis$basis_coefs

train_data <- list(
    input = standard_u_h,
    input_mean = u_h_mean, input_sd = u_h_sd,
    output = fric_arr,
    grounding_line = gl_arr
)

## Save output
if (save_data) {
    setsf <- paste0("sets", sets[1], "-", sets[length(sets)]) # formatC(sets, width=2, flag="0
    dir.create(paste0("training_data/", setsf, "/"))

    saveRDS(train_data, paste0(
        "training_data/", setsf, "/train_data_",
        formatC(first_set, width = 2, flag = "0"), "_", data_date, ".rds"
    ))
}


## Work out how many simulations there are
n_per_set <- dim(thickness_velocity_arr)[1]
n_sims <- n_per_set * length(sets)

## Work out how many observations are needed for a 89/10/1 split
n_train <- floor(n_sims * 0.89)
n_val <- floor(n_sims * 0.1)
n_test <- n_sims - n_train - n_val

last_train_set <- floor(n_train / n_per_set)
leftover_train <- n_train - last_train_set * n_per_set

## Read data from the 2nd set to the last train set
## And save them as training data
for (set in sets[2:last_train_set]) {
    ## Read thickness and velocity data
    thickness_velocity_arr <- readRDS(paste0("./training_data/thickness_velocity_arr_", 
                                            formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
    
    ## Standardise input
    standard_u_h <- (thickness_velocity_arr - u_h_mean) / u_h_sd

    ## Read friction data
    print("Reading friction data...")
    fric_arr <- readRDS(paste0("./training_data/friction_arr_", 
                        formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))

    ## Read grounding line data
    print("Reading grounding line data...")
    gl_arr <- readRDS(paste0("./training_data/gl_arr_", 
                        formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))

    ## Read basis coefficients data
    print("Reading basis coefficients data...")
    fric_basis <- readRDS(paste0("./training_data/friction_basis_", 
                            formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
    basis_coefs <- fric_basis$basis_coefs

    ## Save input-output into a list
    train_data <- list(
        input = standard_u_h,
        input_mean = u_h_mean, input_sd = u_h_sd,
        output = fric_arr,
        grounding_line = gl_arr
    )

    if (save_data) {
        saveRDS(train_data, paste0("./training_data/", setsf, "/train_data_", formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
    }
    
}

## The remaining sets will be split into a bit of training (however much is left to make up 89%), 
## then validation (10%) and test (1%) sets
remain_sets <- setdiff(sets, sets[1:last_train_set])

# for (set in remain_sets) {
inds <- 1:(n_per_set * length(remain_sets))

train_inds <- sample(inds, leftover_train)
val_inds <- sample(setdiff(inds, train_inds), n_val)
test_inds <- setdiff(inds, c(train_inds, val_inds))

## Read and standardise thickness and velocity data
files <- lapply(remain_sets, function(set) paste0("./training_data/thickness_velocity_arr_", 
                                            formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
thickness_velocity_list <- lapply(files, readRDS)
thickness_velocity_arr <- abind(thickness_velocity_list, along = 1)
rm(thickness_velocity_list)

standard_u_h <- (thickness_velocity_arr - u_h_mean) / u_h_sd 

## Read friction data
print("Reading friction data for the remaining sets...")
files <- lapply(remain_sets, function(set) paste0("./training_data/friction_arr_", formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
fric_list <- lapply(files, readRDS)
fric_arr <- abind(fric_list, along = 1)
# rm(fric_list)

## Read grounding line data
print("Reading grounding line data for the remaining sets...")
files <- lapply(remain_sets, function(x) paste0("./training_data/gl_arr_", formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
gl_list <- lapply(files, readRDS)
gl_arr <- abind(gl_list, along = 1)
# rm(gl_list)

## Read basis coefficients data
print("Reading basis coefficients data for the remaining sets...")
files <- lapply(remain_sets, function(x) paste0("./training_data/friction_basis_", formatC(set, width = 2, flag = "0"), "_", data_date, ".rds"))
basis_coefs_list <- lapply(files, function(x) readRDS(x)$basis_coefs)
basis_coefs <- abind(basis_coefs_list, along = 1)
# rm(basis_coefs_list)

train_u_h <- standard_u_h[train_inds, , ,]
val_u_h <- standard_u_h[val_inds, , ,]
test_u_h <- standard_u_h[test_inds, , ,]
# rm(standard_u_h) # to save memory

train_fric <- basis_coefs[train_inds, ]
val_fric <- basis_coefs[val_inds, ]
test_fric <- basis_coefs[test_inds, ]
# rm(basis_coefs) # save memory

train_gl <- as.vector(gl_arr[train_inds, , ,])
val_gl <- as.vector(gl_arr[val_inds, , ,])
test_gl <- as.vector(gl_arr[test_inds, , ,])
# rm(gl_arr) # save memory

## Save input-output into a list
train_data <- list(
    input = train_u_h,
    input_mean = u_h_mean, input_sd = u_h_sd,
    output = fric_arr,
    grounding_line = gl_arr
)

val_data <- list(
input = val_u_h, 
input_mean = u_h_mean, input_sd = u_h_sd,
output = val_fric, 
grounding_line = val_gl
)

test_data <- list(
    input = test_u_h,
    input_mean = u_h_mean, input_sd = u_h_sd,
    output = test_fric,
    grounding_line = test_gl
)

if (save_data) {
    print("Saving output...")
    saveRDS(train_data, file = paste0("./training_data/", setsf, "/train_data_", 
                                        formatC(last_train_set+1, width = 2, flag = "0"), "_", data_date, ".rds"))
    saveRDS(val_data, file = paste0("./training_data/", setsf, "/val_data_", data_date, ".rds"))
    saveRDS(test_data, file = paste0("./training_data/", setsf, "/test_data_", data_date, ".rds"))
}
