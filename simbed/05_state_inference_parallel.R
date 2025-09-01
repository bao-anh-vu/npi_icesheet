### Main file ###

setwd("/home/babv971/SSA_model/CNN/simbed/")

rm(list = ls())

library("parallel")
library("Matrix")
library("qlcMatrix")
library("fastmatrix")
library("expm")
library("R.utils")
library("sp")
library("fields")
library("tidyr")
library("dplyr")
library("ggplot2")
library("matrixStats")
library("mvtnorm")
# library("mvnfast")
library("splines")
library(gridExtra)
library(qs)
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
source("./source/enkf/run_enkf.R")
source("./source/enkf/run_enkf_missing.R")
source("./source/enkf/construct_missing_mat.R")

# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
source("./source/enkf/initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# source("ssa_enkf_plots.R")

use_missing_pattern <- T

## Seed for generating bed
# ssa_seed <- 123
# set.seed(ssa_seed)

run_EnKF <- T
save_enkf_output <- T

## EnKF flags
add_process_noise <- T
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
use_true_velocity <- F # use reference velocity as the initial previous velocity
# use_velocity_dependent_R <- FALSE
# run_analysis <- T
regenerate_ini_ens <- T # generate a fresh initial ensemble
# use_basis_functions <- F
use_true_thickness <- F
use_true_bed <- F
use_true_friction <- F
use_cov_taper <- F # use covariance taper
inflate_cov <- F

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"
Ne <- 500 # Ensemble size
years <- 20 #20 # 40
save_points <- c(1, floor(years/2) + 1, years+1) #c(1, 11, 21)
steps_per_yr <- 52 # 100
n_params <- 1

## Command line args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("At least two arguments must be supplied (sample number and EnKF run number).", call. = FALSE)
}
sample_ind <- as.numeric(args[1]) # sample index from the test set
run <- as.numeric(args[2]) # enkf run number (1-10)

print(paste("Sample index:", sample_ind, "EnKF run number:", run))

## Sample from test set and do state inference
set.seed(2024)
# chosen_test_samples <- sample(1:500, 50)
# set.seed(NULL)
# s <- chosen_test_samples[sample_ind]
s <- sample_ind # test sample

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

## SSA model info
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
# reference <- readRDS(file = paste("./training_data/initial_conds/reference_20220329.rds", sep = ""))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed and friction from NN output
sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

if (use_missing_pattern) {
    data_dir <- paste0("./training_data/", setsf, "/missing")
} else {
    data_dir <- paste0("./training_data/", setsf, "/nonmissing")
}
test_data <- qread(file = paste0(data_dir, "/test_data_", output_date, ".qs"))

true_surface_elev <- test_data$true_surface_elev
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- test_data$true_fric

print("Reading posterior samples...")
if (use_missing_pattern) {
    output_dir <- paste0("./output/posterior/", setsf, "/missing")
} else {
    output_dir <- paste0("./output/posterior/", setsf, "/nonmissing")
}

## Posterior samples
# fric_samples_ls <- readRDS(file = paste0(output_dir, "/fric_post_samples_", output_date, ".rds"))
# bed_samples_ls <- readRDS(file = paste0(output_dir, "/bed_post_samples_", output_date, ".rds"))
# # gl_samples_ls <- readRDS(file = paste0(output_dir, "/gl_post_samples_", output_date, ".rds"))

fric_samples_ls <- qread(file = paste0(output_dir, "/fric_post_samples_", output_date, ".qs"))
bed_samples_ls <- qread(file = paste0(output_dir, "/bed_post_samples_", output_date, ".qs"))
# gl_samples_ls <- readRDS(file = paste0(output_dir, "/gl_post_samples_", output_date, ".rds"))

## Mean prediction
pred_fric <- qread(file = paste0(output_dir, "/pred_fric_", output_date, ".qs"))
pred_bed <- qread(file = paste0(output_dir, "/pred_bed_", output_date, ".qs"))
# pred_gl <- readRDS(pred_gl, file = paste0(output_dir, "/pred_gl_", output_date, ".rds"))

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1 / 3)
# pred_fric <- pred_fric * fric_scale

## Read surface observations
surface_elev <- test_data$input[, , , 1] * test_data$input_sd[1] + test_data$input_mean[1]
velocity <- test_data$input[, , , 2] * test_data$input_sd[2] + test_data$input_mean[2]

if (use_missing_pattern) {
    print("Reading missing patterns...")
    surf_elev_missing_pattern <- readRDS("./training_data/surf_elev_missing_pattern.rds")
    vel_missing_pattern <- readRDS("./training_data/vel_missing_pattern.rds")
    missing_pattern <- list(surface_elev = surf_elev_missing_pattern, vel = vel_missing_pattern)

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

} else {
    missing_pattern <- NULL
}

################################
##      State inference       ##
################################

# for (s in test_samples) {
    
    print(paste("Sample", s))

    if (use_missing_pattern) {
        enkf_output_dir <- paste0("./output/posterior/", setsf, "/missing/sample", sample_ind)
    } else {
        enkf_output_dir <- paste0("./output/posterior/", setsf, "/nonmissing/sample", sample_ind)
    }

    if (!dir.exists(enkf_output_dir)) {
        dir.create(paste0(enkf_output_dir))
    }

    if (use_cov_taper) {
        enkf_output_dir <- paste0(enkf_output_dir, "/taper")
    } else {
        enkf_output_dir <- paste0(enkf_output_dir, "/no_taper")
    }

    if (!dir.exists(enkf_output_dir)) {
        dir.create(paste0(enkf_output_dir))
    } #else { # delete all previously saved plots
    #     unlink(paste0(enkf_output_dir, "/*"))
    # }
    
    fric_s <- pred_fric[, s]
    bed_s <- pred_bed[, s]

    # surface_obs_s <- surface_obs[s,,,]
    surface_elev_s <- surface_elev[s, , ]
    velocity_s <- velocity[s, , ]
    
    if (use_missing_pattern) {
        surface_elev_s[which(missing_pattern$surface_elev == 0)] <- NA
        velocity_s[which(missing_pattern$vel == 0)] <- NA
    }

    # surf_elev_df <- data.frame(x = domain, z = ini_surf_obs) # u = observations$velocity_obs[, 1])
    # polyfit <- loess(z ~ x, data = surf_elev_df, span = 0.05)#$fitted
    # lmfit <- lm(z ~ x, data = surf_elev_df)
    # smoothed_surf_elev <- predict(polyfit, data = surf_elev_df)

    # png("./plots/surf_elev.png")
    # plot(ini_surf_obs, type = "l")
    # dev.off()

    surface_obs_list <- list(surface_elev = surface_elev_s, velocity = velocity_s)

    # 3. Initialise ensemble
    # set.seed(1)
    sample_inds <- sample(1:dim(fric_samples_ls[[s]])[2], n_params, replace = F)
    # set.seed(NULL)
    bed_samples <- bed_samples_ls[[s]][, sample_inds]
    fric_samples <- fric_samples_ls[[s]][, sample_inds]

    if (run_EnKF) {
        ini_ens_list <- list()
        for (p in 1:n_params) {

            if (n_params == 1) {
                bed_sample <- bed_samples
                fric_sample <- fric_samples
            } else {
                bed_sample <- bed_samples[, p]
                fric_sample <- fric_samples[, p]
            }

            print("Initialising ensemble...")
            if (regenerate_ini_ens) {
                ## Bed ##
                if (use_true_bed) {
                    ini_beds <- matrix(rep(true_bed[s, ], Ne), J, Ne)
                } else { ## use samples from parameter posterior
                    ini_beds <- matrix(rep(bed_sample, Ne), J, Ne)
                }

                ## Friction ##
                # print("Simulating friction coefficient...")
                if (use_true_friction) {
                    ini_friction <- matrix(rep((true_fric[s, ]), Ne), J, Ne)
                } else {
                    ini_friction <- matrix(rep(log(fric_sample), Ne), J, Ne)
                    # ini_friction <- log(fric_samples_ls[[s]][, 1:Ne])
                }

                ## Now put the ensemble together
                # print("Calculating ice thicknesses based on simulated beds and observed top surface elevation...")
                
                # for (p in 1:n_params) { # create 20 ensembles from 20 beds

                ## Ice thickness ##
                if (use_true_thickness) {
                    ini_thickness <- matrix(rep(true_thicknesses[s, , 1], Ne), J, Ne)
                    # ini_thickness <- matrix(rep(ssa_steady$current_thickness, Ne), J, Ne)
                } else {

                    ## For the initial surface obs, replace NA values with the last non-NA value
                    ini_surf_obs <- surface_elev_s[, 1]
                    # ini_surf_obs[which(is.na(ini_surf_obs))] <- nonNA[length(nonNA)]
                    
                    # fill NA values with the average of the last 100 grid points #nonNA[length(nonNA)]
                    nonNA <- ini_surf_obs[!is.na(ini_surf_obs)]
                    ini_surf_obs[which(is.na(ini_surf_obs))] <- mean(nonNA[(length(nonNA) - 100):length(nonNA)])


                    ## Process noise parameters
                    ones <- rep(1, length(domain))
                    D <- rdist(domain)
                    l <- 50e3
                    R <- exp_cov(D, l)

                    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
                    L <- t(chol(R))
                    L <- as(L, "dgCMatrix")
                    process_noise_info <- list(corrmat_chol = L, length_scale = l)

                    ini_thickness <- initialise_ice_thickness(domain,
                        n_sims = Ne,
                        surface_obs = ini_surf_obs, # use observed z at t = 0
                        bed = ini_beds[, p],
                        condsim_shelf = F, 
                        process_noise_info = process_noise_info
                    )
                    

                }
                ini_ens <- rbind(
                    ini_thickness,
                    ini_beds,
                    #   matrix(rep(ini_beds[, i], Ne), J, Ne),
                    ini_friction
                )

                ini_ens_list[[p]] <- ini_ens
                # }
                # ini_ens <- rbind(ini_thickness, matrix(rep(ini_beds[, 1], Ne), J, Ne), ini_friction)

                if (save_enkf_output) {
                    qsave(ini_ens_list, file = paste0(enkf_output_dir, "/ini_ens_list_", output_date, ".qs", sep = ""))
                }
            } 
        }
    } else {
        ini_ens_list <- qread(file = paste0(enkf_output_dir, "/ini_ens_list_", output_date, ".qs", sep = ""))
        # ini_ens <- ini_ens[, 1:Ne]
    }

    ## Try plotting one of the ensembles
    # param_ind <- 1
    # ini_state_params <- rbind(ini_ens_list[[param_ind]],
    #                           ini_bed_list[[param_ind]],
    #                           ini_friction_list[[param_ind]])
    # plot_ini_ens(ini_ens, reference, #observations,
    #               print_bed_plot = T, print_friction_plot = T)
    ################################ Initial velocity ##############################

    ini_velocity <- NULL

    if (use_true_velocity) { # use reference velocity as initial velocity
        ini_velocity <- matrix(rep(true_velocities[s, , 1], Ne), J, Ne)
        # ini_velocity <- matrix(rep(ssa_steady$current_velocity, Ne), J, Ne)
    } else { # use a smoothed version of the observed velocity as initial velocity
    
        ini_vel_obs <- velocity_s[, 1]
        velocity_df <- data.frame(x = domain, u = ini_vel_obs) # u = observations$velocity_obs[, 1])
        loessfit <- loess(u ~ x, data = velocity_df, span = 0.05)#$fitted
        
        smoothed_velocity <- predict(loessfit, newdata = data.frame(x = domain))
        
        if (use_missing_pattern) {
            na_ind <- which(is.na(smoothed_velocity))
            polyfit <- lm(u ~ poly(x, 3), data = velocity_df)#$fitted
            smoothed_velocity[na_ind] <- predict(polyfit, newdata = data.frame(x = domain[na_ind]))
        }
        # png(paste0("./plots/velocity_obs.png"), width = 2000, height = 1000, res = 300)
        # plot(velocity_df$u, type = "l")
        # lines(smoothed_velocity, col = "red")
        # dev.off()
        
        ini_velocity <- matrix(rep(smoothed_velocity, Ne), J, Ne)

    }

    ################################################################################
    ##                             EnKF implementation                            ##
    ################################################################################
    if (run_EnKF) {
        ## Process noise parameters
        ones <- rep(1, length(domain))
        D <- rdist(domain)
        l <- 50e3
        R <- exp_cov(D, l)

        # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
        L <- t(chol(R))
        L <- as(L, "dgCMatrix")
        process_noise_info <- list(corrmat_chol = L, length_scale = l)

        enkf_thickness_ls <- list()
        enkf_bed_ls <- list()
        enkf_friction_ls <- list()
        enkf_velocities_ls <- list()

        error_inds <- rep(1, n_params)
        for (p in 1:n_params) {

            ini_ens <- ini_ens_list[[p]]

            test <- try({
                enkf1 <- proc.time()

                enkf_out <- run_enkf_missing_noH(
                domain = domain, years = years,
                steps_per_yr = steps_per_yr,
                ini_thickness = ini_ens[1:J, ],
                ini_bed = ini_ens[J + 1:J, ],
                ini_friction_coef = ini_ens[2 * J + 1:J, ],
                ini_velocity = ini_velocity,
                observations = surface_obs_list,
                missing_pattern = missing_pattern,
                run_analysis = T,
                add_process_noise = add_process_noise,
                use_cov_taper = use_cov_taper,
                inflate_cov = inflate_cov,
                process_noise_info = process_noise_info,
                use_const_measure_error = F
                )

                # enkf_out <- run_enkf_missing_noH(
                # domain = domain, years = years,
                # steps_per_yr = steps_per_yr,
                # ini_thickness = ini_ens[1:J, ],
                # ini_bed = ini_ens[J + 1:J, ],
                # ini_friction_coef = ini_ens[2 * J + 1:J, ],
                # ini_velocity = ini_velocity,
                # observations = surface_obs_list,
                # # missing_pattern = missing_pattern,
                # run_analysis = T,
                # add_process_noise = add_process_noise,
                # use_cov_taper = use_cov_taper,
                # # inflate_cov = inflate_cov,
                # process_noise_info = process_noise_info,
                # use_const_measure_error = F
                # )

            error_inds[p] <- 0 # if run successfully turn off error flag

            enkf_thickness <- enkf_out$ens
            enkf_bed <- enkf_out$bed
            enkf_friction <- enkf_out$friction_coef
            enkf_velocities <- enkf_out$velocities

            # enkf_thickness2 <- enkf_out2$ens
            # enkf_bed2 <- enkf_out2$bed
            # enkf_friction2 <- enkf_out2$friction_coef
            # enkf_velocities2 <- enkf_out2$velocities

            # enkf_thickness3 <- enkf_out3$ens
            # enkf_bed3 <- enkf_out3$bed
            # enkf_friction3 <- enkf_out3$friction_coef
            # enkf_velocities3 <- enkf_out3$velocities

            # K1 <- enkf_out$K
            # K2 <- enkf_out2$K

            # HX1 <- enkf_out$HX
            # HX2 <- enkf_out2$HX

            # yr <- 20
#             # HX1[[yr]][1:5, 1:5] - HX2[[yr]][1:5, 1:5]
#             # K1[[yr]][1:5, 1:5] - K2[[yr]][1:5, 1:5]
            # head(rowMeans(enkf_thickness[[yr]]) - rowMeans(enkf_thickness2[[yr]]))
#             plot(rowMeans(enkf_thickness2[[yr]]) - rowMeans(enkf_thickness[[yr]]))

            enkf_thickness_ls[[p]] <- enkf_thickness
            enkf_bed_ls[[p]] <- enkf_bed
            enkf_friction_ls[[p]] <- enkf_friction
            enkf_velocities_ls[[p]] <- enkf_velocities

            enkf2 <- proc.time()

            }#,
            # error = function(e){error_inds[p] <- 1}
            )
        }

        # print("Running EnKF...")

        # test <- system.time({
        #     enkf_results <- lapply(ini_ens_list, function(ens) {

        #         enkf_out <- run_enkf_missing(
        #         domain = domain, years = years,
        #         steps_per_yr = steps_per_yr,
        #         ini_thickness = ens[1:J, ],
        #         ini_bed = ens[J + 1:J, ],
        #         ini_friction_coef = ens[2 * J + 1:J, ],
        #         ini_velocity = ini_velocity,
        #         observations = surface_obs_list,
        #         missing_pattern = missing_pattern,
        #         run_analysis = T,
        #         add_process_noise = add_process_noise,
        #         use_cov_taper = use_cov_taper,
        #         inflate_cov = inflate_cov,
        #         process_noise_info = process_noise_info
        #         )

        #         enkf_thickness <- enkf_out$ens
        #         enkf_bed <- enkf_out$bed
        #         enkf_friction <- enkf_out$friction_coef
        #         enkf_velocities <- enkf_out$velocities

        #         return(result = list(
        #             thickness = enkf_thickness,
        #             bed = enkf_bed,
        #             friction = enkf_friction,
        #             velocities = enkf_velocities
        #         ))
        # }#,
        # # mc.cores = 2L
        # )
        # })
        
        ## Then concatenate the ensembles
        thickness_concat <- list() # length = number of years
        velocity_concat <- list()

        # for (t in 1:(years+1)) {
        for (ind in 1:length(save_points)) {
            t <- save_points[ind]
            thickness_t <- lapply(enkf_thickness_ls, function(x) x[[t]])
            
            velocity_t <- lapply(enkf_velocities_ls, function(x) x[[t]])
            
            thickness_concat[[ind]] <- do.call(cbind, thickness_t)
            velocity_concat[[ind]] <- do.call(cbind, velocity_t)
        }

        bed_concat <- do.call(cbind, enkf_bed_ls)
        friction_concat <- do.call(cbind, enkf_friction_ls)
    
        if (save_enkf_output) {
            qsave(thickness_concat, file = paste0(enkf_output_dir, "/enkf_thickness_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
            qsave(bed_concat, file = paste0(enkf_output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
            qsave(friction_concat, file = paste0(enkf_output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
            qsave(velocity_concat, file = paste0(enkf_output_dir, "/enkf_velocities_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
            # qsave(error_inds, file = paste0(enkf_output_dir, "/error_inds_sample", s, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
        }
    } else {
        thickness_concat <- qread(file = paste0(enkf_output_dir, "/enkf_thickness_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
        bed_concat <- qread(file = paste0(enkf_output_dir, "/enkf_bed_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
        friction_concat <- qread(file = paste0(enkf_output_dir, "/enkf_friction_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
        velocity_concat <- qread(file = paste0(enkf_output_dir, "/enkf_velocities_sample", sample_ind, "_Ne", Ne, "_", output_date, "_", run, ".qs", sep = ""))
    }

    ################################################################################
    ##                      Plotting the combined ensemble                        ##
    ################################################################################

    # par(mfrow = c(1, 1))
    # combined_enkf_ens <- do.call(cbind, fin_enkf_ens)
    # combined_enkf_beds <- do.call(cbind, fin_enkf_beds)
    # combined_enkf_friction <- do.call(cbind, fin_enkf_friction)
    # combined_enkf_velocities <- do.call(cbind, fin_velocities)

    if (use_missing_pattern) {
        plot_dir <- paste0("./plots/posterior/", setsf, "/missing/sample", sample_ind)
    } else {
        plot_dir <- paste0("./plots/posterior/", setsf, "/nonmissing/sample", sample_ind)
    }

    if (!dir.exists(plot_dir)) {
        dir.create(paste0(plot_dir))
    } 

    if (use_cov_taper) {
        plot_dir <- paste0(plot_dir, "/taper", "_Ne", Ne)
    } else {
        plot_dir <- paste0(plot_dir, "/no_taper", "_Ne", Ne)
    }

    if (!dir.exists(plot_dir)) {
        dir.create(paste0(plot_dir))
    } #else { # delete all previously saved plots
    #     unlink(paste0(plot_dir, "/*"))
    # }

    plot_times <- save_points #seq(1, years + 1, 2)
    ## Ice thickness plot
    if (plot_ice_thickness) {
        print("Plotting ice thickness...")
        # png(paste0(plot_dir, "/thickness_enkf.png"), width = 2000, height = 1000, res = 300)
        # par(mfrow = c(length(plot_times)/2, 2))

        # thickness_plots <- list()
        for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
            t <- plot_times[ind]
            ens_t <- thickness_concat[[ind]]
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
            true_thickness.df <- data.frame(domain = domain / 1000, thickness = true_thicknesses[s, , t])

            # Plot title
            title <- paste("t = ", t - 1, "a")
            thickness_plot <- ggplot() +
                geom_ribbon(
                    data = enkf_thickness.df,
                    aes(domain, ymin = thickness_low, ymax = thickness_up),
                    fill = "red", alpha = 0.2
                ) +
                theme_bw() +
                geom_line(data = true_thickness.df, aes(domain, thickness)) +
                geom_line(data = enkf_thickness.df, aes(domain, mean_thickness), colour = "red") +
                xlab("Domain (km)") +
                ylab("Ice thickness (m)") +
                ggtitle(title) +
                theme(plot.title = element_text(
                    hjust = 0.95, vjust = 0.5, face = "bold",
                    margin = margin(t = 20, b = -30)
                ))

            # thickness_plots[[ind]] <- thickness_plot
            png(paste0(plot_dir, "/thickness_enkf_", t - 1, ".png"), width = 2000, height = 1000, res = 300)
            print(thickness_plot)
            # grid.arrange(grobs = thickness_plots, ncol = 2)
            dev.off()
        }


        # grid.arrange(grobs = thickness_plots, ncol = 2)
        # dev.off()
    }


    ## Velocity plot

    if (plot_velocity) {
        print("Plotting velocity...")
        # vel_plots <- list()
        for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
            t <- plot_times[ind]
            ens_t <- velocity_concat[[ind]]
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
                lq = thickness_low, uq = thickness_up
            )
            true_velocities.df <- data.frame(domain = domain / 1000, velocity = true_velocities[s, , t])
            # Plot title
            title <- paste("t = ", t - 1, "a")
            velocity_plot <- ggplot() +
                geom_ribbon(
                    data = enkf_velocities.df,
                    aes(domain, ymin = vel_low, ymax = vel_up),
                    fill = "red", alpha = 0.2
                ) +
                theme_bw() +
                geom_line(data = true_velocities.df, aes(domain, velocity)) +
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

            png(paste0(plot_dir, "/velocity_enkf", t - 1, ".png"), width = 2000, height = 1000, res = 300)
            print(velocity_plot)
            # grid.arrange(grobs = vel_plots, ncol = 2)
            dev.off()
            # print(thickness_plot)
        }
    }


    # plot_times <- seq(1, years + 1, 5)
    ## Bed plot
    if (plot_bed) {
        print("Plotting bed...")
        png(paste0(plot_dir, "/bed.png"), width = 2000, height = 1000, res = 300)
        # matplot(domain/1000, combined_bg_beds, type = "l", col = "lightgrey",
        #         xlab = "Domain (km)", ylab = "Elevation (m)")
        
        # enkf_velocities.df <- data.frame(
        #         domain = domain / 1000, mean_velocity = rowMeans(ens_t),
        #         lq = thickness_low, uq = thickness_up
        #     )
        
        matplot(domain / 1000, bed_concat,
            type = "l", col = "lightpink",
            ylab = "Bed elevation (m)", xlab = "Domain (km)"
        )
        # matlines(domain/1000, combined_pf_beds, col = "skyblue2")
        # lines(domain/1000, reference$bedrock, col = "black", lwd = 1.5)
        lines(domain / 1000, true_bed[s, ], col = "black", lwd = 1.5)
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
        print("Plotting friction...")
        png(paste0(plot_dir, "/friction.png"), width = 2000, height = 1000, res = 300)
        plot_region <- 1:J # 500:1000 #(2*J+1):(3*J)
        # matplot(domain[plot_region]/1000, 10^combined_bg_friction[plot_region, ], type = "l", col = "lightgrey",
        #         xlab = "Domain (km)", ylab = "Friction coefficient (unit)")
        matplot(domain[plot_region] / 1000, exp(friction_concat[plot_region, ]),
            type = "l", col = "lightpink",
            ylab = "Friction (unit)", xlab = "Domain (km)"
        )
        # matlines(domain[plot_region]/1000, 10^combined_pf_friction[plot_region, ], col = "skyblue2")
        # lines(domain[plot_region]/1000, reference$friction_coef[plot_region], col = "black", lwd = 1.5)
        lines(domain / 1000, exp(true_fric[s, ]), col = "black", lwd = 1.5)
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

    ## RMSE (time series)
    rmse <- function(estimated, true) {
        stopifnot(length(estimated) == length(true))
        sum(sqrt((estimated - true)^2))
    }

    enkf.rmse <- c()
    for (i in 1:length(save_points)) {
        t <- save_points[i]
        enkf.rmse[i] <- rmse(rowMeans(thickness_concat[[i]]), true_thicknesses[s, , t])
    }
    plot(save_points, enkf.rmse,
        type = "o", col = "red", xlab = "Time (years)", ylab = "RMSE",
        main = paste0("RMSE of ice thickness over time for sample", s)
    )

# }
