### Main file ###

setwd("/home/babv971/SSA_model/CNN/real_data/")

rm(list = ls())

library(parallel)
library(Matrix)
library(qlcMatrix)
library(fastmatrix)
library(expm)
library(R.utils)
library(sp)
library(fields)
library(tidyr)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(mvtnorm)
# library("mvnfast")
library(splines)
library(gridExtra)
library(qs)
library(abind)
# library("fda")

source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/surface_elev.R")
source("./source/enkf/get_obs.R")
source("./source/enkf/initialise_ens.R")
source("./source/enkf/ssa_plot_ini_ens.R")
source("./source/enkf/propagate.R")
source("./source/enkf/obs_operator.R")
# source("./source/enkf/run_enkf.R")
source("./source/enkf/run_enkf_missing.R")
source("./source/enkf/construct_missing_mat.R")
source("./source/simulate_friction.R")

# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
source("./source/enkf/initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# source("ssa_enkf_plots.R")

## Seed for generating bed
# ssa_seed <- 123
# set.seed(ssa_seed)

run_EnKF <- T
save_enkf_output <- T
# save_bg_output <- F

## EnKF flags
add_process_noise <- T
# use_true_velocity <- F # use reference velocity as the initial previous velocity
regenerate_ini_ens <- T # generate a fresh initial ensemble
use_basis_functions <- F
# use_true_thickness <- F
# use_true_bed <- F
# use_true_friction <- F
use_cov_taper <- F # use covariance taper
inflate_cov <- F
correct_model_discrepancy <- T
correct_velocity_discrepancy <- T # can correct velocity discrepancy as well, if F will just correct surface elevation discrepancy
avg_over_time <- T
# use_posterior_samples <- F

## Presets
data_date <- "20241111" # "20230518"
output_date <- "20241111" # "20240518"
Ne <- 10 # Ensemble size
years <- 11 # 40
steps_per_yr <- 100
n_params <- 1 # number of beds

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

## Command line args
args <- commandArgs(trailingOnly = TRUE)
sample_type <- as.character(args[1]) # "prior" or "posterior"
p <- as.numeric(args[2]) # sample index from the test set

## Read bed and friction from NN output
sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

output_dir <- paste0("./output/cnn/", setsf, "/")

## Directories
if (correct_model_discrepancy) {
    if (avg_over_time) {
        pred_output_dir <- paste0(output_dir, "pred/discr_avg/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr_avg/")
    } else {
        pred_output_dir <- paste0(output_dir, "pred/discr/")
        plot_dir <- paste0("./plots/cnn/", setsf, "/pred/discr/")
    }
} else {
    pred_output_dir <- paste0(output_dir, "pred/")
    plot_dir <- paste0("./plots/cnn/", setsf, "/pred/")
}

# Load ice sheet at steady state
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Load real data
surf_elev_data <- qread(file = "./data/surface_elev/surf_elev_mat.qs")
velocity_data <- qread(file = "./data/velocity/vel_smoothed.qs")

## "Mask" locations beyond the GL in the velocity data
se_na_locs <- which(is.na(surf_elev_data[, 1]))
velocity_data[se_na_locs, ] <- NA

## Mask velocity for 2012 as well because it's way too noisy and upsets the filter
velocity_data[, 1:10] <- NA

## Also "mask" the last year for leave-one-out validation
surf_elev_data[, ncol(surf_elev_data)] <- NA
velocity_data[, ncol(velocity_data)] <- NA

vel_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/vel_discr_smooth_", data_date, ".qs"))
se_discr_mat <- qread(file = paste0("./data/discrepancy/", setsf, "/se_discr_smooth_", data_date, ".qs"))

# Extrapolate discrepancies to the 11th year (discrepancy is constant in time so just repeat last year)
# se_discr_mat <- cbind(se_discr_mat, se_discr_mat[, ncol(se_discr_mat)]) # repeat last column to make 11 years
# vel_discr_mat <- cbind(vel_discr_mat, vel_discr_mat[, ncol(vel_discr_mat)]) # repeat last column to make 11 years

model_discrepancy <- list(
    surface_elev = se_discr_mat,
    velocity = vel_discr_mat
)

## Load measurement noise info
vel_err_sd <- qread(file = paste0("./data/velocity/vel_err_sd_", data_date, ".qs"))
se_err_sd <- rep(0.5, length(domain)) #qsave(file = paste0("./data/surface_elev/se_err_sd_", data_date, ".qs"))
measurement_noise_info <- list(
    vel_err_sd = vel_err_sd,
    se_err_sd = se_err_sd
)

# Load physical parameters (SMB, basal melt rate, etc.)
params <- qread(file = paste0("./data/training_data/", "/phys_params_", data_date, ".qs"))
# params$ab <- 0

# test_data <- qread(file = paste0(data_dir, "/test_data_", data_date, ".qs"))

# true_surface_elev <- test_data$true_surface_elev
# true_thicknesses <- test_data$true_thickness_test
# true_velocities <- test_data$true_velocity_test

# true_bed <- test_data$true_bed
# true_fric <- test_data$true_fric

if (sample_type == "posterior") {
    # sample_type <- "posterior"
    # Posterior samples
    print("Reading posterior samples...")
    fric_samples_ls <- qread(file = paste0(pred_output_dir, "fric_samples_real_", data_date, ".qs"))
    bed_samples_ls <- qread(file = paste0(pred_output_dir, "bed_samples_real_", data_date, ".qs"))
} else if (sample_type == "prior") {
    # Prior samples
    # print("Simulating from prior...")
    # set.seed(2025)

    # n_samples <- 1000
    # bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
    # L <- t(chol(bed_prior$cov))
    # u_mat <- matrix(rnorm(nrow(L) * n_samples), nrow = nrow(L), ncol = n_samples)
    # mean_mat <- matrix(rep(bed_prior$mean, n_samples), nrow = nrow(bed_prior$mean), ncol = n_samples)
    # bed_samples_ls <- mean_mat + L %*% u_mat

    # fric_samples_ls <- simulate_friction2(
    #     nsim = n_samples, domain = domain
    # )

    # qsave(bed_samples_ls, file = paste0(pred_output_dir, "prior_bed_samples_", data_date, ".qs"))
    # qsave(fric_samples_ls, file = paste0(pred_output_dir, "prior_fric_samples_", data_date, ".qs"))

    print("Reading prior samples...")
    bed_samples_ls <- qread(file = paste0(pred_output_dir, "prior_bed_samples_basis_", data_date, ".qs"))
    fric_samples_ls <- qread(file = paste0(pred_output_dir, "prior_fric_samples_basis_", data_date, ".qs"))

} else {
    stop("sample_type must be either 'prior' or 'posterior'")
}



# if (use_missing_pattern) {
#     print("Reading missing patterns...")
#     surf_elev_missing_pattern <- qread("./data/surface_elev/missing_pattern.qs")
#     vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")
#     missing_pattern <- list(surface_elev = surf_elev_missing_pattern, vel = vel_missing_pattern)

# vel_mp_c <- c(vel_missing_pattern)
# vel_mp_df <- data.frame(gridpt = rep(1:J, ncol(vel_missing_pattern)),
#                         nonmissing = vel_mp_c,
#                         year = rep(0:(ncol(vel_missing_pattern)-1), each = J))

# ## Plot missing velocity pattern over space-time
# vel_st_plot <- ggplot(vel_mp_df) +
#     geom_tile(aes(x = gridpt, y = year, fill = factor(nonmissing))) +
#     # scale_fill_distiller(palette = "BuPu", direction = 1) +
#     theme_bw() +
#     labs(x = "Grid Point", y = "Year", fill = "Non-missing") +
#     ggtitle("Thwaites Glacier Velocity Over Time") + theme(plot.title = element_text(hjust = 0.5))


# png("./plots/temp/vel_mp.png")
# print(vel_st_plot)
# dev.off()
# }
# } else {
#     missing_pattern <- NULL
# }


################################
##      State inference       ##
################################

## Sample from test set and do state inference
# n_test_samples <- dim(test_data$input)[1]
# test_samples <- 1 # seq(100, 500, 100) # sample(1:n_test_samples, 1) # sample index

# s <- test_samples[1]

# if (use_missing_pattern) {
enkf_output_dir <- paste0(pred_output_dir, "enkf/")
# } else {
#     enkf_output_dir <- paste0("./output/cnn/", setsf, "/nonmissing/sample", s)
# }

# if (!dir.exists(enkf_output_dir)) {
#     dir.create(paste0(enkf_output_dir))
# }

# if (use_cov_taper) {
#     enkf_output_dir <- paste0(enkf_output_dir, "/taper")
# } else {
#     enkf_output_dir <- paste0(enkf_output_dir, "/no_taper")
# }

if (!dir.exists(enkf_output_dir)) {
    dir.create(paste0(enkf_output_dir))
} # else { # delete all previously saved plots
#     unlink(paste0(output_dir, "/*"))
# }t

# for (s in test_samples) {

## Replace NA values with the last non-NA value
# ini_surf_obs <- surface_elev_s[, 1]
# nonNA <- ini_surf_obs[!is.na(ini_surf_obs)]
# ini_surf_obs[which(is.na(ini_surf_obs))] <- nonNA[length(nonNA)]

surface_obs_list <- list(surface_elev = surf_elev_data, velocity = velocity_data)

# 3. Initialise ensemble
# bed_samples <- bed_samples_ls[, 1:n_params]
# fric_samples <- fric_samples_ls[, 1:n_params]

if (run_EnKF) {
    bed_sample <- bed_samples_ls[, p]
    fric_sample <- fric_samples_ls[, p]

    # cat("Bed:", head(bed_sample), "\n")
    # cat("Fric:", head(fric_sample))

    print("Initialising ensemble...")
    # if (regenerate_ini_ens) {
    ## Bed ##
    ini_beds <- matrix(rep(bed_sample, Ne), J, Ne)

    ## Friction ##
    ini_friction <- matrix(rep(log(fric_sample), Ne), J, Ne)

    # ## Process noise parameters
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 50e3
    R <- exp_cov(D, l)

    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    process_noise_info <- list(corrmat_chol = L, length_scale = l)

    ini_thickness_mean <- ssa_steady$current_thickness #+ rnorm(J, 0, 10)
    h_mat <- matrix(rep(ini_thickness_mean, Ne), nrow = length(ini_thickness_mean), ncol = Ne)

    h_sd <- c(rep(5, length(domain)))

    h_sd_mat <- matrix(rep(h_sd, Ne), nrow = length(h_sd), ncol = Ne)
    Zmat <- matrix(rnorm(length(domain) * Ne), nrow = length(domain), ncol = Ne)
    h_noise <- h_sd_mat * (L %*% Zmat)
    ini_thickness <- h_mat + h_noise


    ini_ens <- rbind(
        ini_thickness,
        ini_beds,
        ini_friction
    )
}

## Plot initial ice thickness ensembles
png(paste0(plot_dir, "/enkf/ini_thickness_ens.png"), width = 2000, height = 3000, res = 300)
par(mfrow = c(3, 1))

matplot(ini_thickness, type = "l")

matplot(ini_beds, type = "l")

matplot(exp(ini_friction), type = "l")
dev.off()

################################ Initial velocity ##############################

ini_velocity <- matrix(rep(ssa_steady$current_velocity, Ne), J, Ne)

################################################################################
##                             EnKF implementation                            ##
################################################################################
if (run_EnKF) {

    enkf_thickness_ls <- list()
    enkf_bed_ls <- list()
    enkf_friction_ls <- list()
    enkf_velocities_ls <- list()

    # error_inds <- rep(1, n_params)
    # for (p in 1:n_params) {

    #     ini_ens <- ini_ens_list[[p]]

    #     test <- try({
    #         enkf1 <- proc.time()

    #         enkf_out <- run_enkf_missing_noH(
    #         domain = domain, years = years,
    #         steps_per_yr = steps_per_yr,
    #         ini_thickness = ini_ens[1:J, ],
    #         ini_bed = ini_ens[J + 1:J, ],
    #         ini_friction_coef = ini_ens[2 * J + 1:J, ],
    #         ini_velocity = ini_velocity,
    #         observations = surface_obs_list,
    #         missing_pattern = missing_pattern,
    #         # run_analysis = F,
    #         add_process_noise = add_process_noise,
    #         use_cov_taper = use_cov_taper,
    #         inflate_cov = inflate_cov,
    #         process_noise_info = process_noise_info
    #         )

    #         # enkf_out2 <- run_enkf_missing(
    #         # domain = domain, years = years,
    #         # steps_per_yr = steps_per_yr,
    #         # ini_thickness = ini_ens[1:J, ],
    #         # ini_bed = ini_ens[J + 1:J, ],
    #         # ini_friction_coef = ini_ens[2 * J + 1:J, ],
    #         # ini_velocity = ini_velocity,
    #         # observations = surface_obs_list,
    #         # missing_pattern = missing_pattern,
    #         # # run_analysis = F,
    #         # add_process_noise = add_process_noise,
    #         # use_cov_taper = use_cov_taper,
    #         # inflate_cov = inflate_cov,
    #         # process_noise_info = process_noise_info
    #         # )

    #     error_inds[p] <- 0 # if run successfully turn off error flag

    #     enkf_thickness <- enkf_out$ens
    #     enkf_bed <- enkf_out$bed
    #     enkf_friction <- enkf_out$friction_coef
    #     enkf_velocities <- enkf_out$velocities

    #     # enkf_thickness2 <- enkf_out2$ens
    #     # enkf_bed2 <- enkf_out2$bed
    #     # enkf_friction2 <- enkf_out2$friction_coef
    #     # enkf_velocities2 <- enkf_out2$velocities

    #     # png("./plots/temp/ens_thickness.png")
    #     #  plot(rowMeans(enkf_thickness[[1]]), type = "l")
    #     #  lines(rowMeans(enkf_thickness2[[1]]), col = "red")
    #     # dev.off()

    #     enkf_thickness_ls[[p]] <- enkf_thickness
    #     enkf_bed_ls[[p]] <- enkf_bed
    #     enkf_friction_ls[[p]] <- enkf_friction
    #     enkf_velocities_ls[[p]] <- enkf_velocities

    #     enkf2 <- proc.time()

    #     }#,
    #     # error = function(e){error_inds[p] <- 1}
    #     )
    # }

    print("Running EnKF...")

    test <- system.time({
        # enkf_results <- lapply(ini_ens_list, function(ens) {

        enkf_out <- run_enkf_missing_noH(
            domain = domain, years = years,
            steps_per_yr = steps_per_yr,
            ini_thickness = ini_ens[1:J, ],
            ini_bed = ini_ens[J + 1:J, ],
            ini_friction_coef = ini_ens[2 * J + 1:J, ],
            ini_velocity = ini_velocity,
            observations = surface_obs_list,
            phys_params = params,
            missing_pattern = missing_pattern,
            correct_model_discrepancy = T,
            model_discrepancy = model_discrepancy,
            run_analysis = T,
            add_process_noise = F,
            use_cov_taper = use_cov_taper,
            inflate_cov = inflate_cov,
            process_noise_info = process_noise_info,
            measurement_noise_info = measurement_noise_info
        )

        enkf_thickness_ls <- enkf_out$ens
        enkf_se_ls <- enkf_out$surf_elev_ens
        enkf_bed_ls <- enkf_out$bed
        enkf_friction_ls <- enkf_out$friction_coef
        enkf_velocities_ls <- enkf_out$velocities

        # return(result = list(
        #     thickness = enkf_thickness,
        #     bed = enkf_bed,
        #     friction = enkf_friction,
        #     velocities = enkf_velocities
        # ))
        # },
        # mc.cores = 2L
        # )
    })

    # enkf_thickness_ls <- lapply(enkf_results, function(r) r$thickness)
    # enkf_bed_ls <- lapply(enkf_results, function(r) r$bed)
    # enkf_friction_ls <- lapply(enkf_results, function(r) r$friction)
    # enkf_velocities_ls <- lapply(enkf_results, function(r) r$velocities)

    ## Then concatenate the ensembles

    # thickness_concat <- list() # length = number of years
    # velocity_concat <- list()

    # for (t in 1:(years+1)) {
    #     thickness_t <- lapply(enkf_thickness_ls, function(x) x[[t]])

    #     velocity_t <- lapply(enkf_velocities_ls, function(x) x[[t]])
    #     thickness_concat[[t]] <- do.call(cbind, thickness_t)

    #     velocity_concat[[t]] <- do.call(cbind, velocity_t)
    # }

    # bed_concat <- do.call(cbind, enkf_bed_ls)
    # friction_concat <- do.call(cbind, enkf_friction_ls)

    thickness_concat <- enkf_thickness_ls
    se_concat <- enkf_se_ls
    bed_concat <- enkf_bed_ls
    friction_concat <- enkf_friction_ls
    velocity_concat <- enkf_velocities_ls

    if (save_enkf_output) {
        qsave(thickness_concat, file = paste0(enkf_output_dir, "/enkf_thickness_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
        qsave(se_concat, file = paste0(enkf_output_dir, "/enkf_surf_elev_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
        qsave(bed_concat, file = paste0(enkf_output_dir, "/enkf_bed_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
        qsave(friction_concat, file = paste0(enkf_output_dir, "/enkf_friction_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
        qsave(velocity_concat, file = paste0(enkf_output_dir, "/enkf_velocities_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
        # qsave(error_inds, file = paste0(enkf_output_dir, "/error_inds_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
    }
} else {
    thickness_concat <- qread(file = paste0(enkf_output_dir, "/enkf_thickness_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
    bed_concat <- qread(file = paste0(enkf_output_dir, "/enkf_bed_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
    friction_concat <- qread(file = paste0(enkf_output_dir, "/enkf_friction_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
    velocity_concat <- qread(file = paste0(enkf_output_dir, "/enkf_velocities_", sample_type, "_p", p, "_", output_date, ".qs", sep = ""))
}

################################################################################
##                      Plotting the combined ensemble                        ##
################################################################################

# par(mfrow = c(1, 1))
# combined_enkf_ens <- do.call(cbind, fin_enkf_ens)
# combined_enkf_beds <- do.call(cbind, fin_enkf_beds)
# combined_enkf_friction <- do.call(cbind, fin_enkf_friction)
# combined_enkf_velocities <- do.call(cbind, fin_velocities)

plot_times <- c(1, 10) + 1 # +1 because time index starts from 1, and we want to exclude the initial condition
## Ice thickness plot
if (plot_ice_thickness) {
    # png(paste0(plot_dir, "/thickness_enkf.png"), width = 2000, height = 1000, res = 300)
    # par(mfrow = c(length(plot_times)/2, 2))

    # thickness_plots <- list()
    # for (ind in 1:length(plot_times)) {
    for (t in plot_times) {
        # t <- plot_times[ind]
        ens_t <- thickness_concat[[t]]
        state_var <- 1 / (ncol(ens_t) - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
        thickness_low <- rowMeans(ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(state_var)[1:J])
        thickness_up <- rowMeans(ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(state_var)[1:J])

        # Grounding line
        # ref.GL <- reference$grounding_line[t]
        # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
        # GL <- GL / J * domain[length(domain)] / 1000

        enkf_thickness.df <- data.frame(
            domain = domain / 1000, mean_thickness = rowMeans(ens_t),
            lq = thickness_low, uq = thickness_up
        )
        # true_thickness.df <- data.frame(domain = domain / 1000, thickness = true_thicknesses[s, , t])

        # Plot title
        title <- paste("t = ", t - 1, "a")
        thickness_plot <- ggplot() +
            geom_ribbon(
                data = enkf_thickness.df,
                aes(domain, ymin = thickness_low, ymax = thickness_up),
                fill = "red", alpha = 0.2
            ) +
            theme_bw() +
            # geom_line(data = true_thickness.df, aes(domain, thickness)) +
            geom_line(data = enkf_thickness.df, aes(domain, mean_thickness), colour = "red") +
            xlab("Domain (km)") +
            ylab("Ice thickness (m)") +
            ggtitle(title) +
            theme(plot.title = element_text(
                hjust = 0.95, vjust = 0.5, face = "bold",
                margin = margin(t = 20, b = -30)
            ))

        # thickness_plots[[ind]] <- thickness_plot
        png(paste0(plot_dir, "/enkf/thickness_enkf_", t - 1, ".png"), width = 2000, height = 1000, res = 300)
        print(thickness_plot)
        # grid.arrange(grobs = thickness_plots, ncol = 2)
        dev.off()
    }


    # grid.arrange(grobs = thickness_plots, ncol = 2)
    # dev.off()
}


## Velocity plot

if (plot_velocity) {
    # vel_plots <- list()
    # for (ind in 1:length(plot_times)) {
    for (t in plot_times) {
        # t <- plot_times[ind]
        ens_t <- velocity_concat[[t]]
        state_var <- 1 / (ncol(ens_t) - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
        vel_low <- rowMeans(ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(state_var)[1:J])
        vel_up <- rowMeans(ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(state_var)[1:J])

        # Grounding line
        # ref.GL <- reference$grounding_line[t]
        # GL <- gl_migrate(H = enkf_means[1:J, t], b = rowMeans(ens[(J+1):(2*J), ]))
        # GL <- GL / J * domain[length(domain)] / 1000
        # bg.GL <- gl_migrate(H = bg_state_means[1:J, t], b = bg_state_means[(J+1):(2*J), t])
        # bg.GL <- bg.GL / J * domain[length(domain)] / 1000
        #

        enkf_velocities.df <- data.frame(
            domain = domain / 1000, mean_velocity = rowMeans(ens_t),
            lq = vel_low, uq = vel_up
        )
        # true_velocities.df <- data.frame(domain = domain / 1000, velocity = true_velocities[s, , t])
        # Plot title
        title <- paste("t = ", t - 1, "a")
        velocity_plot <- ggplot() +
            geom_ribbon(
                data = enkf_velocities.df,
                aes(domain, ymin = vel_low, ymax = vel_up),
                fill = "red", alpha = 0.25
            ) +
            theme_bw() +
            # geom_line(data = true_velocities.df, aes(domain, velocity)) +
            geom_line(data = enkf_velocities.df, aes(domain, mean_velocity), colour = "red") +
            # geom_line(data = enkf_thickness.df, aes(domain, lq), colour = "red") +
            # geom_line(data = enkf_thickness.df, aes(domain, uq), colour = "blue") +

            # geom_line(data = bg_thickness.df,
            #           aes(domain, bg_thickness.df[, t]), colour = "gray") +
            # geom_line(data = reference_thickness.df,
            #           aes(domain, reference_thickness.df[, t]), colour = "black") +
            # geom_line(data = enkf_thickness.df,
            #           aes(domain, enkf_thickness.df[, t]), colour = "red") +
            # geom_vline(xintercept = ref.GL[ind], linetype="dashed",
            #            color = "black", size = 0.5) +
            # geom_vline(xintercept = bg.GL[ind], linetype="dashed",
            #            color = "gray", size = 0.5) +
            # geom_vline(xintercept = enkf.GL[ind], linetype="dashed",
            #            color = "red", size = 0.5) +
            # coord_cartesian(xlim = c(0, 450), ylim = c(200, 3000)) +
            xlab("Domain (km)") +
            ylab("Velocity (m/a)") +
            ggtitle(title) +
            theme(plot.title = element_text(
                hjust = 0.95, vjust = 0.5, face = "bold",
                # margin = margin(t = 20, b = -30)
                margin = margin(t = 0, b = 0)
            ))

        png(paste0(plot_dir, "/enkf/velocity_enkf", t - 1, ".png"), width = 2000, height = 1000, res = 300)
        print(velocity_plot)
        # grid.arrange(grobs = vel_plots, ncol = 2)
        dev.off()
        # print(thickness_plot)
    }
}

## Bed plot
if (plot_bed) {
    png(paste0(plot_dir, "/enkf/bed.png"), width = 2000, height = 1000, res = 300)
    # matplot(domain/1000, combined_bg_beds, type = "l", col = "lightgrey",
    #         xlab = "Domain (km)", ylab = "Elevation (m)")

    enkf_velocities.df <- data.frame(
        domain = domain / 1000, mean_velocity = rowMeans(ens_t),
        lq = thickness_low, uq = thickness_up
    )

    matplot(domain / 1000, bed_concat,
        type = "l", col = "lightpink",
        ylab = "Bed elevation (m)", xlab = "Domain (km)"
    )
    # matlines(domain/1000, combined_pf_beds, col = "skyblue2")
    # lines(domain/1000, reference$bedrock, col = "black", lwd = 1.5)
    # lines(domain / 1000, true_bed[s, ], col = "black", lwd = 1.5)
    lines(domain / 1000, rowMeans(bed_concat), col = "red", lwd = 1.5)
    # lines(domain/1000, rowMeans(combined_pf_beds), col = "royalblue", lwd = 1.5)
    legend("bottomleft",
        legend = c("True bed", "Predicted bed"), # "Background ensemble",
        # "EnKF ensemble", "EnKF mean"),
        col = c(
            "black", # "grey", "pink",
            "red"
        ), lty = 1
    )
    dev.off()
}

# Friction plot
if (plot_friction) {
    png(paste0(plot_dir, "/enkf/friction.png"), width = 2000, height = 1000, res = 300)
    plot_region <- 1:J # 500:1000 #(2*J+1):(3*J)
    # matplot(domain[plot_region]/1000, 10^combined_bg_friction[plot_region, ], type = "l", col = "lightgrey",
    #         xlab = "Domain (km)", ylab = "Friction coefficient (unit)")
    matplot(domain[plot_region] / 1000, exp(friction_concat[plot_region, ]),
        type = "l", col = "lightpink",
        ylab = "Friction (unit)", xlab = "Domain (km)"
    )
    # matlines(domain[plot_region]/1000, 10^combined_pf_friction[plot_region, ], col = "skyblue2")
    # lines(domain[plot_region]/1000, reference$friction_coef[plot_region], col = "black", lwd = 1.5)
    # lines(domain / 1000, exp(true_fric[s, ]), col = "black", lwd = 1.5)
    lines(domain[plot_region] / 1000, rowMeans(exp(friction_concat[plot_region, ])), col = "red", lwd = 1.5)
    # lines(domain[plot_region]/1000, rowMeans(10^combined_pf_friction[plot_region, ]), col = "royalblue", lwd = 1.5)
    legend("topleft",
        legend = c("True friction", "Predicted friction"),
        # "Background ensemble",
        # "PF ensemble", "PF mean",
        # "EnKF ensemble", "EnKF mean"),
        col = c(
            "black", # "lightgrey",
            # "skyblue2", "royalblue",
            #  "lightpink",
            "red"
        ), lty = 1, cex = 0.7
    )
    dev.off()
}

# ## RMSE (time series)
# rmse <- function(estimated, true) {
#     stopifnot(length(estimated) == length(true))
#     sum(sqrt((estimated - true)^2))
# }

# enkf.rmse <- c()
# for (t in 1:years) {
#     enkf.rmse[t] <- rmse(rowMeans(thickness_concat[[t + 1]]), true_thicknesses[s, , t + 1])
# }
# plot(1:years, enkf.rmse,
#     type = "o", col = "red", xlab = "Time (years)", ylab = "RMSE",
#     main = paste0("RMSE of ice thickness over time for sample", s)
# )

# }
