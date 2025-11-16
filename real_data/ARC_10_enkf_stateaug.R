### EnKF with state augmentation (Gillet-Chaulet 2020)

setwd("~/SSA_model/CNN/real_data/")

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
library(FRK)
library(qs)
library(mgcv)
# library("fda")

source("./source/sim_params.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/fit_basis.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/surface_elev.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/enkf/get_obs.R")
source("./source/enkf/initialise_ens.R")
source("./source/enkf/ssa_plot_ini_ens.R")
source("./source/enkf/propagate.R")
source("./source/enkf/obs_operator.R")
# source("./source/enkf/run_enkf_stateaug.R")
source("./source/enkf/run_enkf_stateaug_missing.R")
source("./source/enkf/construct_missing_mat.R")
source("./source/enkf/initialise_ice_thickness.R")

## Flags
run_EnKF <- T
save_enkf_output <- T
use_missing_pattern <- T
use_basal_melt_data <- F

## EnKF flags
add_process_noise <- F
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
# use_true_velocity <- F # use reference velocity as the initial previous velocity
# use_velocity_dependent_R <- FALSE
# run_analysis <- T
# regenerate_ini_ens <- T # generate a fresh initial ensemble
use_basis_funs <- T
use_cov_taper <- F # use covariance taper
inflate_cov <- F

## Presets
data_date <- "20241111" # "20230518"
output_date <- "20241111" # "20240518"
Ne <- 10#0 # 0 # Ensemble size
years <- 10 # 40
# save_points <- c(1, 5, 10) #c(1, floor(years / 2) + 1, years + 1)
steps_per_yr <- 100

output_dir <- paste0("./output/aug_enkf/")
# data_dir <- paste0("./data/training_data/", setsf, "/missing")
plot_dir <- paste0("./plots/aug_enkf/")

if (!dir.exists(plot_dir)) { dir.create(paste0(plot_dir)) }

if (!dir.exists(output_dir)) { dir.create(paste0(output_dir))} 
# else { # delete all previously saved plots
#     unlink(paste0(output_dir, "/*"))
# }

## Plot settings
plot_ice_thickness <- T
plot_velocity <- T
plot_bed <- T
plot_friction <- T

## SSA model info
ssa_steady <- qread(file = paste0("./data/training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Read bed and friction from NN output
# sets <- 51:100 # 10
# setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Surface elevation
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only
vel_mat <- qread("./data/velocity/vel_smoothed.qs")

print("Reading missing patterns...")
surf_elev_missing_pattern <- qread("./data/surface_elev/missing_pattern.qs")
vel_missing_pattern <- qread("./data/velocity/missing_pattern.qs")
# missing_patterns <- abind(list(surf_elev_missing_pattern, vel_missing_pattern), along = 3)

## Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
bed_obs_chosen <- bed_obs_df %>% filter(chosen == 1)

################################
##      State inference       ##
################################

## Physical params
params <- list(
  secpera = 31556926, # seconds per annum
  n = 3.0, # exponent in Glen's flow law
  rho_i = 917.0, # ice density
  rho_w = 1028.0, # sea water density
  g = 9.81 # gravity constant
  # A = 4.227e-25, #1.4579e-25, # flow rate parameter
)

params$m <- 1 / params$n
params$B <- 0.55 * 1e6 * params$secpera^params$m
params$A <- params$B^(-params$n)

## SMB data
smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)
params$as <- smb_avg # surface accumulation rate (m/s)

if (use_basal_melt_data) {
  melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
  # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
  avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
  melt_nonmissing <- which(!is.na(avg_melt_rate))
  avg_melt_rate[1:(melt_nonmissing[1]-1)] <- -1 #seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
  avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) #mean(avg_melt_rate[melt_nonmissing])
  avg_melt_rate <- - avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
} else {
  avg_melt_rate <- rep(0, J)
}
params$ab <- avg_melt_rate # melt rate (m/s)

# Generate samples of bed and friction
print("Simulating friction coefficients...")
fric_sims <- simulate_friction2(
    nsim = Ne, domain = domain
)

## Simulate beds
print("Simulating beds...")
bed_sim_output <- simulate_bed(Ne,
    domain = domain,
    obs_locations = bed_obs_chosen$ind,
    obs = bed_obs_chosen$bed_elev
)

bed_sims <- bed_sim_output$sims
bed_mean <- bed_sim_output$mean

if (use_basis_funs) {
    n_fric_basis <- 120
    n_bed_basis <- 150

    ## Fit basis to log(friction)
    friction_basis <- fit_friction_basis(
        nbasis = n_fric_basis,
        domain = domain,
        fric_arr = t(fric_sims),
        log_transform = T,
        lengthscale = 4e3
    )

    ## Fit basis to bedrock

    bed_basis <- fit_bed_basis(
        nbasis = n_bed_basis, domain = domain,
        bed_arr = t(bed_sims), lengthscale = 5e3
    )
    # bed_fit <- list(mean = bed_mean, basis = bed_basis)

    friction_ens <- t(friction_basis$fitted_values)
    bed_ens <- t(bed_basis$fitted_values)
} else {
    friction_ens <- fric_sims
    bed_ens <- bed_sims
}

## Initialise ice thickness ensemble
# surface_elev_s <- surface_elev[s, , ]
# velocity_s <- velocity[s, , ]

# if (use_missing_pattern) {
#     surface_elev_s[missing_pattern$surface_elev == 0] <- NA
#     velocity_s[missing_pattern$vel == 0] <- NA
#     surface_elev_s <- matrix(surface_elev_s, nrow = dim(surface_elev)[2], ncol = dim(surface_elev)[3])
#     velocity_s <- matrix(velocity_s, nrow = dim(velocity)[2], ncol = dim(velocity)[3])
# }

# 3. Initialise ensemble

print("Calculating ice thicknesses based on simulated beds and observed top surface elevation...")
## Ice thickness ##
se_grounded <- na.omit(rowMeans(surf_elev_mat)) # Use surface elevation in the year 2000 to initialise ice thickness
ini_thickness <- matrix(0, nrow = J, ncol = Ne)
for (i in 1:Ne) {
    H_ini <- se_grounded - bed_ens[1:length(se_grounded), i]
    length_shelf <- J - length(se_grounded)
    # H_gl <- - bed_sim[(gl_ind+1)] * params$rho_w / params$rho_i #- 100 # minus an offset to satisfy flotation condition
    H_shelf <- rep(500, length_shelf) # seq(from = H_gl, to = 500, length.out = length_shelf)
    ini_thickness[, i] <- c(H_ini, H_shelf)
}
# ini_thickness <- matrix(rep(ssa_steady$current_thickness, Ne), J, Ne)

print("Initialising ensemble...")
ini_ens <- rbind(ini_thickness, bed_ens, friction_ens)

## Try plotting one of the ensembles
# param_ind <- 1
# ini_state_params <- rbind(ini_ens_list[[param_ind]],
#                           ini_bed_list[[param_ind]],
#                           ini_friction_list[[param_ind]])
# plot_ini_ens(ini_ens, reference, #observations,
#               print_bed_plot = T, print_friction_plot = T)

################################ Initial velocity ##############################

years <- ncol(vel_mat)
vel_curr <- rowMeans(vel_mat[, (years-5):years]) # average over the last 10 years

## Smooth the velocity out with gam()
vel_df <- data.frame(domain = domain, vel = vel_curr)
vel_data <- vel_df %>% filter(!is.na(vel))
vel_missing <- vel_df %>% filter(is.na(vel))

gam_fit <- gam(vel ~ s(domain, k = 20), data = vel_data)

# Predict at missing locations (works even if outside observed range)
vel_pred <- predict(
  gam_fit,
  newdata = data.frame(domain = vel_missing$domain)
)
vel_missing$vel <- vel_pred

vel_curr_smooth <- rbind(vel_data, vel_missing) %>%
  arrange(domain) %>%
  pull(vel)

ini_velocity <- matrix(rep(vel_curr_smooth, Ne), J, Ne)

# ini_velocity <- matrix(rep(ssa_steady$current_velocity, Ne), J, Ne)
################################################################################
##                             EnKF implementation                            ##
################################################################################
if (run_EnKF) {
    # ## Process noise parameters
    # ones <- rep(1, length(domain))
    # D <- rdist(domain)
    # l <- 50e3
    # R <- exp_cov(D, l)

    # # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    # L <- t(chol(R))
    # L <- as(L, "dgCMatrix")
    # process_noise_info <- list(corrmat_chol = L, length_scale = l)

    # enkf_thickness_ls <- list()
    # enkf_bed_ls <- list()
    # enkf_friction_ls <- list()
    # enkf_velocities_ls <- list()
    # ini_ens <- ini_ens_list[[i]]

# test <- solve_ssa_nl(
#             domain = domain,
#             bedrock = bed_ens[,2],
#             friction_coef = exp(friction_ens[,2]) * 1e6 * params$secpera^(1 / params$n),
#             phys_params = params,
#             ini_velocity = ssa_steady$current_velocity,#ini_velocity[, 1],
#             ini_thickness = ssa_steady$current_thickness, #ini_ens[1:J, 1],
#             years = years, 
#             steps_per_yr = 100,
#             add_process_noise = F
#         )

    missing_pattern <- list(surface_elev = surf_elev_missing_pattern, vel = vel_missing_pattern)

    ## Apply missing/masking pattern to observations
    surf_elev_mat[surf_elev_missing_pattern == 0] <- NA
    vel_mat[vel_missing_pattern == 0] <- NA
    surface_obs_list <- list(surface_elev = surf_elev_mat, velocity = vel_mat)


    enkf1 <- proc.time()

    enkf_out <- run_enkf_missing(
        domain = domain, years = years,
        steps_per_yr = steps_per_yr,
        phys_params = params,
        # ini_thickness = ini_ens[1:J, ],
        # ini_bed = ini_ens[J + 1:J, ],
        # ini_friction_coef = ini_ens[2 * J + 1:J, ],
        ini_ens = ini_ens,
        ini_velocity = ini_velocity,
        observations = surface_obs_list,
        missing_pattern = missing_pattern,
        run_analysis = T,
        add_process_noise = add_process_noise,
        use_cov_taper = use_cov_taper,
        # inflate_cov = inflate_cov,
        process_noise_info = process_noise_info,
        use_const_measure_error = T
    )

    

    # enkf_out <- run_enkf(
    #     domain = domain, years = years,
    #     steps_per_yr = steps_per_yr,
    #     ini_thickness = ini_ens[1:J, ],
    #     ini_bed = ini_ens[J + 1:J, ],
    #     ini_friction_coef = ini_ens[2 * J + 1:J, ],
    #     ini_velocity = ini_velocity,
    #     observations = surface_obs_list,
    #     # missing_pattern = missing_pattern,
    #     run_analysis = T,
    #     add_process_noise = add_process_noise,
    #     use_cov_taper = use_cov_taper,
    #     process_noise_info = process_noise_info,
    #     use_const_measure_error = F
    # )

    enkf_thickness <- lapply(enkf_out$ens, function(x) x[1:J, ])
    enkf_bed <- lapply(enkf_out$ens, function(x) x[J + 1:J, ])
    enkf_friction <- lapply(enkf_out$ens, function(x) x[2 * J + 1:J, ])
    enkf_velocities <- enkf_out$velocities

    # enkf_thickness <- enkf_thickness[save_points]
    # enkf_bed <- enkf_bed[save_points]
    # enkf_friction <- enkf_friction[save_points]
    # enkf_velocities <- enkf_velocities[save_points]

    enkf2 <- proc.time()

    if (save_enkf_output) {
        qsave(enkf_thickness, file = paste0(output_dir, "/enkf_thickness_", output_date, ".qs", sep = ""))
        qsave(enkf_bed, file = paste0(output_dir, "/enkf_bed_", output_date, ".qs", sep = ""))
        qsave(enkf_friction, file = paste0(output_dir, "/enkf_friction_", output_date, ".qs", sep = ""))
        qsave(enkf_velocities, file = paste0(output_dir, "/enkf_velocities_", output_date, ".qs", sep = ""))
    }
} else {
    enkf_thickness <- qread(file = paste0(output_dir, "/enkf_thickness_", output_date, ".qs", sep = ""))
    enkf_bed <- qread(file = paste0(output_dir, "/enkf_bed_", output_date, ".qs", sep = ""))
    enkf_friction <- qread(file = paste0(output_dir, "/enkf_friction_", output_date, ".qs", sep = ""))
    enkf_velocities <- qread(file = paste0(output_dir, "/enkf_velocities_", output_date, ".qs", sep = ""))
}

################################################################################
##                      Plotting the combined ensemble                        ##
################################################################################

# par(mfrow = c(1, 1))
# combined_enkf_ens <- do.call(cbind, fin_enkf_ens)
# combined_enkf_beds <- do.call(cbind, fin_enkf_beds)
# combined_enkf_friction <- do.call(cbind, fin_enkf_friction)
# combined_enkf_velocities <- do.call(cbind, fin_velocities)


if (use_cov_taper) {
    plot_dir <- paste0(plot_dir, "/taper")
} else {
    plot_dir <- paste0(plot_dir, "/no_taper")
}

if (!dir.exists(plot_dir)) {
    dir.create(paste0(plot_dir))
} # else { # delete all previously saved plots
#     unlink(paste0(plot_dir, "/*"))
# }

plot_times <- seq(1, years + 1, 2)

## Ice thickness plot
if (plot_ice_thickness) {
    print("Plotting ice thickness...")
    # png(paste0(plot_dir, "/thickness_enkf.png"), width = 2000, height = 1000, res = 300)
    # par(mfrow = c(length(plot_times)/2, 2))

    # thickness_plots <- list()
    for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
        t <- plot_times[ind]
        ens_t <- enkf_thickness[[ind]]
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
    print("Plotting velocity...")
    # vel_plots <- list()
    for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
        t <- plot_times[ind]
        ens_t <- enkf_velocities[[ind]]
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
    }
}
# plot_times <- seq(1, years + 1, 5)

## Bed plot
if (plot_bed) {
    print("Plotting bed...")
    for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
        t <- plot_times[ind]
        ens_t <- enkf_bed[[ind]]
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
    print("Plotting friction...")
    for (ind in 1:length(plot_times)) {
        # for (t in plot_times) {
        t <- plot_times[ind]
        ens_t <- enkf_friction[[ind]]
        state_var <- as.matrix(1 / (Ne - 1) * tcrossprod(ens_t - rowMeans(ens_t))) # diag(enkf_covmats[[t]])

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
# }
