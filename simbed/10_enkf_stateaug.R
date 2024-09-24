### EnKF with state augmentation (Gillet-Chaulet 2020)

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
library(FRK)
# library("fda")

source("./source/sim_params.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/fit_basis.R")
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/enkf/get_obs.R")
source("./source/enkf/initialise_ens.R")
source("./source/enkf/ssa_plot_ini_ens.R")
source("./source/enkf/propagate.R")
source("./source/enkf/obs_operator.R")
source("./source/enkf/run_enkf_stateaug.R")
source("./source/enkf/surface_elev.R")

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
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
use_true_velocity <- F # use reference velocity as the initial previous velocity
# use_velocity_dependent_R <- FALSE
# run_analysis <- T
regenerate_ini_ens <- T # generate a fresh initial ensemble
use_basis_funs <- T
use_true_thickness <- F
use_true_bed <- F
use_true_friction <- F
use_cov_taper <- T # use covariance taper
log_transform <- T

## Presets
data_date <- "20220329" # "20230518"
output_date <- "20240320" # "20240518"
Ne <- 100#0 # Ensemble size
years <- 20 # 40
steps_per_yr <- 52 # 100
# n_params <- 1 # 20 #number of beds
# n_bed_obs <- 100
# smoothing_factor <- 0.5 # for the PF

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

## SSA model info
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))
reference <- readRDS(file = paste("./training_data/initial_conds/reference_20220329.rds", sep = ""))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed and friction from NN output
sets <- 1:50 # 10
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

data_dir <- paste0("./training_data/", setsf)

if (use_basis_funs) {
    output_dir <- paste0("./output/stateaug/", setsf, "/basis")    
} else {
    output_dir <- paste0("./output/stateaug/", setsf, "/no_basis")
}

# if (!dir.exists(output_dir)) {
#     dir.create(paste0(output_dir))
# } else { # delete previously created output dir
#     unlink(paste0(output_dir, "/*"))
# }

## Read test data
test_data <- readRDS(file = paste0(data_dir, "/test_data_", output_date, ".rds"))

true_surface_elevs <- test_data$true_surface_elevs_test
true_thicknesses <- test_data$true_thickness_test
true_velocities <- test_data$true_velocity_test

true_bed <- test_data$true_bed
true_fric <- test_data$true_fric

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1 / 3)
# pred_fric <- pred_fric * fric_scale

## Read surface observations
surface_elev <- test_data$input[, , , 1] * test_data$input_sd[1] + test_data$input_mean[1]
velocity <- test_data$input[, , , 2] * test_data$input_sd[2] + test_data$input_mean[2]

## Read bed observations
bed_obs <- readRDS(file = paste0("./training_data/bed_obs_", output_date, ".rds"))

################################
##      State inference       ##
################################

## Sample from test set and do state inference
n_test_samples <- dim(test_data$input)[1]
test_samples <- 1 # seq(100, 500, 100) # sample(1:n_test_samples, 1) # sample index

# NEED TO GENERATE ENSEMBLES OF BED AND FRICTION HERE
sim_param_list <- sim_params(
      nsims = Ne, domain = ssa_steady$domain,
      bed_obs = bed_obs
    )

if (use_basis_funs) {
    ## Fit basis to log(friction)
    nbasis <- 150

    friction_basis <- fit_friction_basis(
        nbasis = nbasis,
        domain = domain,
        fric_arr = sim_param_list$friction,
        log_transform = log_transform
    )

    ## Fit basis to bedrock

    bed_arr <- sim_param_list$bedrock
    bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, bed_arr = bed_arr)
    # bed_fit <- list(mean = bed_mean, basis = bed_basis)

    friction_ens <- t(friction_basis$fitted_values)
    bed_ens <- t(bed_basis$fitted_values)

} else {
    friction_ens <- log(t(sim_param_list$friction))
    bed_ens <- t(sim_param_list$bedrock)
}

s <- 1
# for (s in test_samples) {
    print(paste("Sample", s))

    # surface_obs_s <- surface_obs[s,,,]
    surface_elev_s <- surface_elev[s, , ]
    velocity_s <- velocity[s, , ]

    surface_obs_list <- list(surface_elev = surface_elev_s, velocity = velocity_s)

    # 3. Initialise ensemble

    if (run_EnKF) {
        ini_ens_list <- list()

        print("Initialising ensemble...")
        if (regenerate_ini_ens) {
            ## Bed ##
            if (use_true_bed) {
                ini_beds <- matrix(rep(true_bed[s, ], Ne), J, Ne)
            } else { ## use samples from parameter posterior
                ini_beds <- bed_ens
            }

            ## Friction ##
            print("Simulating friction coefficient...")
            if (use_true_friction) {
                ini_friction <- matrix(rep((true_fric[s, ]), Ne), J, Ne)
            } else {
                # ini_friction <- matrix(rep(log(fric_sample), Ne), J, Ne)
                ini_friction <- friction_ens
            }

            ## Now put the ensemble together
            print("Calculating ice thicknesses based on simulated beds and observed top surface elevation...")
            
            # for (p in 1:n_params) { # create 20 ensembles from 20 beds

            ## Ice thickness ##
            if (use_true_thickness) {
                ini_thickness <- matrix(rep(true_thicknesses[s, , 1], Ne), J, Ne)
                # ini_thickness <- matrix(rep(ssa_steady$current_thickness, Ne), J, Ne)
            } else {
                # ini_thickness <- matrix(0, J, Ne)
                # for (i in 1:Ne) { ## parallelise
                #     ini_thickness[, i] <- initialise_ice_thickness(domain,
                #     n_sims = 1,
                #     surface_obs = surface_elev_s[, 1], # use observed z at t = 0
                #     bed = ini_beds[, i],
                #     condsim_shelf = F
                #     )
                # }

                test <- system.time({
                    
                        bed_list <- lapply(1:ncol(bed_ens), function(i) bed_ens[, i])

                        ini_thickness_ls <- lapply(bed_list, initialise_ice_thickness, 
                        domain = domain,
                        n_sims = 1,
                        surface_obs = surface_elev_s[, 1], # use observed z at t = 0
                        # mc.cores = 10L
                        )
                    
                    ini_thickness <- matrix(unlist(ini_thickness_ls), nrow = J, ncol = Ne)
                    })
            }
            ini_ens <- rbind(
                ini_thickness,
                ini_beds,
                #   matrix(rep(ini_beds[, i], Ne), J, Ne),
                ini_friction
            )

            ini_ens_list[[s]] <- ini_ens
            # }
            # ini_ens <- rbind(ini_thickness, matrix(rep(ini_beds[, 1], Ne), J, Ne), ini_friction)

            if (save_enkf_output) {
                saveRDS(ini_ens_list, file = paste0(output_dir, "/ini_ens_list_", output_date, "_Ne", Ne, ".rds", sep = ""))
            }
        } else {
            ini_ens_list <- readRDS(file = paste0(output_dir, "/ini_ens_list_", output_date, "_Ne", Ne, ".rds", sep = ""))
            # ini_ens <- ini_ens[, 1:Ne]
        }
    # }

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
        velocity_df <- data.frame(x = domain, u = velocity_s[, 1]) # u = observations$velocity_obs[, 1])
        smoothed_velocity <- loess(u ~ x, data = velocity_df, span = 0.05)$fitted
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


        ini_ens <- ini_ens_list[[s]]
        enkf1 <- proc.time()

        enkf_out <- run_enkf(
            domain = domain, years = years,
            steps_per_yr = steps_per_yr,
            ini_thickness = ini_ens[1:J, ],
            ini_bed = ini_ens[J + 1:J, ],
            ini_friction_coef = ini_ens[2 * J + 1:J, ],
            ini_velocity = ini_velocity,
            observations = surface_obs_list,
            run_analysis = T,
            add_process_noise = add_process_noise,
            use_cov_taper = use_cov_taper,
            process_noise_info = process_noise_info,
        )

        enkf_thickness <- lapply(enkf_out$ens, function(x) x[1:J, ])
        enkf_bed <- lapply(enkf_out$ens, function(x) x[J + 1:J, ])
        enkf_friction <- lapply(enkf_out$ens, function(x) x[2 * J + 1:J, ])
        # enkf_bed <- enkf_out$bed
        # enkf_friction <- enkf_out$friction_coef
        enkf_velocities <- enkf_out$velocities

        enkf2 <- proc.time()

        if (save_enkf_output) {
            saveRDS(enkf_thickness, file = paste0(output_dir, "/enkf_thickness_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
            saveRDS(enkf_bed, file = paste0(output_dir, "/enkf_bed_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
            saveRDS(enkf_friction, file = paste0(output_dir, "/enkf_friction_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
            saveRDS(enkf_velocities, file = paste0(output_dir, "/enkf_velocities_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
        }
    } else {
        enkf_thickness <- readRDS(file = paste0(output_dir, "/enkf_thickness_sample_", "_Ne", Ne, output_date, ".rds", sep = ""))
        enkf_bed <- readRDS(file = paste0(output_dir, "/enkf_bed_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
        enkf_friction <- readRDS(file = paste0(output_dir, "/enkf_friction_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
        enkf_velocities <- readRDS(file = paste0(output_dir, "/enkf_velocities_sample_", output_date, "_Ne", Ne, ".rds", sep = ""))
    }

    ################################################################################
    ##                      Plotting the combined ensemble                        ##
    ################################################################################

    # par(mfrow = c(1, 1))
    # combined_enkf_ens <- do.call(cbind, fin_enkf_ens)
    # combined_enkf_beds <- do.call(cbind, fin_enkf_beds)
    # combined_enkf_friction <- do.call(cbind, fin_enkf_friction)
    # combined_enkf_velocities <- do.call(cbind, fin_velocities)

    print("Saving plots...")
    if (use_basis_funs) {
        plot_dir <- paste0("./plots/stateaug/", setsf, "/basis")
    } else {
        plot_dir <- paste0("./plots/stateaug/", setsf, "/no_basis")
    }

    if (!dir.exists(plot_dir)) {
        dir.create(paste0(plot_dir))
    } else { # delete all previously saved plots
        unlink(paste0(plot_dir, "/*"))
    }

    plot_times <- seq(1, years + 1, 2)
    ## Ice thickness plot
    if (plot_ice_thickness) {
        # png(paste0(plot_dir, "/thickness_enkf.png"), width = 2000, height = 1000, res = 300)
        # par(mfrow = c(length(plot_times)/2, 2))

        # thickness_plots <- list()
        # for (ind in 1:length(plot_times)) {
        for (t in plot_times) {
            # t <- plot_times[ind]
            ens_t <- enkf_thickness[[t]]
            state_var <- 1 / (Ne - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
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
        # vel_plots <- list()
        # for (ind in 1:length(plot_times)) {
        for (t in plot_times) {
            # t <- plot_times[ind]
            ens_t <- enkf_velocities[[t]]
            state_var <- 1 / (Ne - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
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


        # png(paste0(plot_dir, "/velocity.png"), width = 2000, height = 1000, res = 300)
        # # matplot(domain/1000, combined_bg_velocities[1:J, ], type = "l", col = "lightgrey",
        # #         xlab = "Domain (km)", ylab = "Velocity (m/a)")
        # # matlines(domain/1000, combined_pf_velocities[1:J, ], lty = 2, col = "skyblue2")
        # matplot(domain/1000, enkf_velocities[[years+1]], lty = 2, type = "l", col = "lightpink",
        #   ylab = "Velocity (m/a)", xlab = "Domain (km)", main = paste0("Velocity at year ", years))
        # # lines(domain/1000, rowMeans(combined_bg_velocities[1:J, ]), col = "grey")
        # # lines(domain/1000, reference$all_velocities[, years+1], col = "black", lwd = 1.5)
        # lines(domain/1000, true_velocities[s, , years+1], col = "black", lwd = 1.5)
        # lines(domain/1000, rowMeans(enkf_velocities[[years+1]]), col = "red", lwd = 1.5)
        # # lines(domain/1000, rowMeans(combined_pf_velocities[1:J, ]), col = "royalblue", lwd = 1.5)
        # legend("bottomright", legend = c("True thickness",
        #                                   # "Background ensemble",
        #                                  "EnKF ensemble", "EnKF mean"), #, "PF ensemble", "PF mean"),
        #        col = c("black", #"grey",
        #        "lightpink", "red"), #, "skyblue2", "royalblue"),
        #        lty = 1, cex = 0.7)
        # dev.off()
    }

    ## Bed plot
    if (plot_bed) {

        for (t in plot_times) {
            # t <- plot_times[ind]
            ens_t <- enkf_bed[[t]]
            state_var <- 1 / (Ne - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
            lower <- rowMeans(ens_t[1:J, ]) + qnorm(0.025) * sqrt(diag(state_var)[1:J])
            upper <- rowMeans(ens_t[1:J, ]) + qnorm(0.975) * sqrt(diag(state_var)[1:J])

            enkf_bed.df <- data.frame(
                domain = domain / 1000, mean = rowMeans(ens_t),
                lq = lower, uq = upper
            )
            true_bed.df <- data.frame(domain = domain / 1000, val = true_bed[s, ])
            # Plot title
            title <- paste("t = ", t - 1, "a")
            bed_plot <- ggplot() +
                geom_ribbon(
                    data = enkf_bed.df,
                    aes(domain, ymin = lower, ymax = upper),
                    fill = "red", alpha = 0.2
                ) +
                theme_bw() +
                geom_line(data = true_bed.df, aes(domain, val)) +
                geom_line(data = enkf_bed.df, aes(domain, mean), colour = "red") +
                xlab("Domain (km)") +
                ylab("Bed elevation (m)") +
                ggtitle(title) +
                theme(plot.title = element_text(
                    hjust = 0.95, vjust = 0.5, face = "bold",
                    # margin = margin(t = 20, b = -30)
                    margin = margin(t = 0, b = 0)
                ))

            png(paste0(plot_dir, "/bed_enkf", t - 1, ".png"), width = 2000, height = 1000, res = 300)
            print(bed_plot)
            # grid.arrange(grobs = vel_plots, ncol = 2)
            dev.off()
        }
    }

    # Friction plot
    if (plot_friction) {

        for (t in plot_times) {
            # t <- plot_times[ind]
            ens_t <- enkf_friction[[t]]
            state_var <- 1 / (Ne - 1) * tcrossprod(ens_t - rowMeans(ens_t)) # diag(enkf_covmats[[t]])
            
            ## Need to sample from Gaussian posterior of the log(friction) here,
            ## then exponentiate to get the friction values
            ## and lastly find quantiles

            log_fric_post_samples <- rmvnorm(1000, rowMeans(ens_t), state_var)
            fric_post_samples <- exp(log_fric_post_samples)
            fric_post_mean <- colMeans(fric_post_samples)
            lower <- apply(fric_post_samples, 2, quantile, probs = 0.05)
            upper <- apply(fric_post_samples, 2, quantile, probs = 0.95)

            enkf_friction.df <- data.frame(
                domain = domain / 1000, mean = fric_post_mean,
                lq = lower, uq = upper
            )
            true_friction.df <- data.frame(domain = domain / 1000, val = exp(true_fric[s, ]))
            # Plot title
            title <- paste("t = ", t - 1, "a")
            friction_plot <- ggplot() +
                geom_ribbon(
                    data = enkf_friction.df,
                    aes(domain, ymin = lower, ymax = upper),
                    fill = "red", alpha = 0.2
                ) +
                theme_bw() +
                geom_line(data = true_friction.df, aes(domain, val)) +
                geom_line(data = enkf_friction.df, aes(domain, mean), colour = "red") +
                xlab("Domain (km)") +
                ylab("Friction (unit)") +
                ggtitle(title) +
                theme(plot.title = element_text(
                    hjust = 0.95, vjust = 0.5, face = "bold",
                    # margin = margin(t = 20, b = -30)
                    margin = margin(t = 0, b = 0)
                ))

            png(paste0(plot_dir, "/friction_enkf", t - 1, ".png"), width = 2000, height = 1000, res = 300)
            print(friction_plot)
            # grid.arrange(grobs = vel_plots, ncol = 2)
            dev.off()
        }
    }

    
}
