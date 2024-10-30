## Assemble EnKF-CNN output

library(abind)
library(qs)

setwd("/home/babv971/SSA_model/CNN/simbed/")

sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

set.seed(2024)
chosen_test_samples <- sample(1:500, 50) # choose 50 samples from the test set
set.seed(NULL)

s <- 7
Ne <- 500
years <- 20
save_points <- c(1, floor(years/2) + 1, years+1) #c(1, 11, 21)
output_date <- "20240320" # "20240518"
use_missing_pattern <- T
# add_process_noise <- T
use_cov_taper <- F # use covariance taper
inflate_cov <- F

runs <- 1:10

if (use_missing_pattern) {
    enkf_output_dir <- paste0("./output/posterior/", setsf, "/missing/sample", s)
} else {
    enkf_output_dir <- paste0("./output/posterior/", setsf, "/nonmissing/sample", s)
}

if (use_cov_taper) {
    enkf_output_dir <- paste0(enkf_output_dir, "/taper")
} else {
    enkf_output_dir <- paste0(enkf_output_dir, "/no_taper")
}

thickness_files <- lapply(runs, function(r) paste0(enkf_output_dir, "/enkf_thickness_sample", s, "_Ne", Ne, "_", output_date, "_", r, ".qs", sep = ""))
bed_files <- lapply(runs, function(r) paste0(enkf_output_dir, "/enkf_bed_sample", s, "_Ne", Ne, "_", output_date, "_", r, ".qs", sep = ""))
fric_files <- lapply(runs, function(r) paste0(enkf_output_dir, "/enkf_friction_sample", s, "_Ne", Ne, "_", output_date, "_", r, ".qs", sep = ""))
velocity_files <- lapply(runs, function(r) paste0(enkf_output_dir, "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date, "_", r, ".qs", sep = ""))

thickness_concat <- list()
bed_concat <- list()
fric_concat <- list()
velocity_concat <- list()

for (p in 1:length(save_points)) {
    thickness_ens <- lapply(thickness_files, qread) 
    thickness_year <- lapply(thickness_ens, function(x) as.matrix(x[[p]]))
    thickness_concat[[p]] <- abind(thickness_year, along = 2)

    velocity_ens <- lapply(velocity_files, qread)
    velocity_year <- lapply(velocity_ens, function(x) as.matrix(x[[p]]))
    velocity_concat[[p]] <- abind(velocity_year, along = 2)
}

bed_ens <- lapply(bed_files, qread)
bed_concat <- abind(bed_ens, along = 2)

fric_ens <- lapply(fric_files, qread)
fric_concat <- abind(fric_ens, along = 2)

print("Saving output...")
qsave(thickness_concat, file = paste0(enkf_output_dir, "/enkf_thickness_sample", s, "_Ne", Ne, "_", output_date, ".qs", sep = ""))
qsave(bed_concat, file = paste0(enkf_output_dir, "/enkf_bed_sample", s, "_Ne", Ne, "_", output_date, ".qs", sep = ""))
qsave(fric_concat, file = paste0(enkf_output_dir, "/enkf_friction_sample", s, "_Ne", Ne, "_", output_date, ".qs", sep = ""))
qsave(velocity_concat, file = paste0(enkf_output_dir, "/enkf_velocities_sample", s, "_Ne", Ne, "_", output_date, ".qs", sep = ""))    

## Remove individual files
rm_thickness_files <- lapply(thickness_files, file.remove)
rm_bed_files <- lapply(bed_files, file.remove)
rm_fric_files <- lapply(fric_files, file.remove)
rm_velocity_files <- lapply(velocity_files, file.remove)

