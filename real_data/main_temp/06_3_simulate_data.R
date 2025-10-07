## Generate parameters from the prior

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library(parallel)
library(Matrix)
library(R.utils)
library(fields)
library(dplyr)
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(abind)
library(ggplot2)
library(gridExtra)
library(qs)
# library(mgcv)

# source("./source/sim_params.R")
# source("./source/run_sims.R")
source("./source/sim_obs.R")
source("./source/cond_sim_gp.R")
source("./source/get_ini_thickness.R")
source("./source/process_sim_results.R")
# source("./source/fit_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
# source("./source/create_ref.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/get_surface_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

## Some flags
# regenerate_sims <- T
# save_sims <- T
# log_transform <- T
use_basal_melt_data <- T
constrain_gl <- F
use_basis <- F

## Directory for training data
train_data_dir <- "./data/training_data"

## Presets
data_date <- "20241111"
# N <- 10 # 00 # number of simulations per set
# set <- 1 #commandArgs(trailingOnly = TRUE)
sets <- 11:20 #50 #:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])
warmup <- 0
years <- 11 # number of years data is collected (not including initial condition)

#########################
## Set physical params ##
#########################

params <- list(
    secpera = 31556926, # seconds per annum
    n = 3.0, # exponent in Glen's flow law
    rho_i = 917.0, # ice density
    rho_w = 1028.0, # sea water density
    g = 9.81 # gravity constant
    # A = 4.227e-25, #1.4579e-25, # flow rate parameter
)

params$m <- 1 / params$n
params$B <- 0.6 * 1e6 * params$secpera^params$m
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
    avg_melt_rate[1:(melt_nonmissing[1] - 1)] <- -1 # seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
    avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) # mean(avg_melt_rate[melt_nonmissing])
    avg_melt_rate <- -avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
} else {
    avg_melt_rate <- rep(0, J)
}
params$ab <- avg_melt_rate # melt rate (m/s)

# Save physical params
qsave(params, file = paste0(train_data_dir, "/phys_params_", data_date, ".qs"))

##########################
##      Read data       ##
##########################

# Ice sheet at steady state (to use as initial conditions)
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

# Surface elevation data
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

# Bed observations
# bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))
# bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]
bed_obs_df <- qread(file = paste0("./data/bedmap/bed_obs_df_all.qs"))
bed_obs_chosen <- bed_obs_df

# if (constrain_gl) {
#   ## Impose an exta condition at GL: bed = surface elevation / (rho_w/rho_i + 1)

#   steady_se <- ssa_steady$current_top_surface
#   steady_gl <- ssa_steady$grounding_line[length(ssa_steady$grounding_line)]
#   gl_ind <- sum(domain/1e3 < steady_gl)
#   bed_gl <- steady_se[gl_ind] / (1 - params$rho_w / params$rho_i)
#   bed_gl_df <- data.frame(ind = gl_ind, loc = domain[gl_ind], bed_elev = bed_gl, bed_sd = 0, interval = NA, chosen = 1)
#   bed_obs_df <- rbind(bed_obs_df, bed_gl_df) %>% arrange(ind)
# }
# bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]
# bed_obs_val <- bed_obs_df[bed_obs_df$chosen == 0, ]

## Bedmachine data to compare
bedmachine_data <- qread(paste0("./data/bedmachine/flowline_bedmachine.qs"))
bedmachine <- bedmachine_data$bed_avg

## Observed grounding line position
gl_pos <- qread(file = paste0("./data/grounding_line/gl_pos.qs"))
gl_ind <- gl_pos$ind

##########################################################
##  Generate observations of surface elev and velocity  ##
##########################################################

flags <- c()

# bed_sim_list <- list()
# fric_sim_list <- list()

t1 <- proc.time()

for (i in 1:length(sets)) {
    ## Read fitted parameters
    set <- sets[i]
    setf <- formatC(set, width = 2, flag = "0")
    cat("Generating simulations for set", set, "\n")



    if (use_basis) {
        bed_basis <- qread(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
        friction_basis <- qread(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
        bed_mean <- bed_basis$mean
        N <- nrow(bed_basis$fitted_values) # number of simulations in this set
        ## Pair the bed-friction observations into a list
        param_list <- lapply(1:N, function(r) {
            list(
                friction = exp(friction_basis$fitted_values[r, ]),
                bedrock = bed_basis$fitted_values[r, ] + bed_mean
            )
        })
    } else {
        fric_sims <- qread(file = paste0(train_data_dir, "/fric_sims_", setf, "_", data_date, ".qs"))
        bed_sims <- qread(file = paste0(train_data_dir, "/bed_sims_", setf, "_", data_date, ".qs"))

        N <- ncol(bed_sims) # number of simulations in this set

        param_list <- lapply(1:N, function(r) {
            list(
                friction = fric_sims[, r], #* 1e6 * params$secpera^(1 / params$n),
                bedrock = bed_sims[, r] #+ bed_mean
            )
        })
    }

    test <- try(
        sim_results <- sim_obs(
            param_list = param_list,
            domain = domain,
            phys_params = params,
            years = years,
            warmup = warmup,
            ini_thickness = ssa_steady$current_thickness,
            ini_velocity = ssa_steady$current_velocity,
            smb = smb_avg,
            basal_melt = avg_melt_rate
        )
    )

    ## Dispose of any simulations that failed
    bad_sims <- sim_results$bad_sims
    if (length(bad_sims) > 0) {
        cat("Some simulations failed in set", set, "\n")
        good_sims <- sim_results$results[-bad_sims]
    } else {
        good_sims <- sim_results$results
    }

    ## Extract observations & corresponding parameters from the simulations
    generated_data <- process_sim_results(sims = good_sims)

    surface_obs_arr_s <- generated_data$surface_obs_arr
    friction_arr_s <- generated_data$friction_arr
    gl_arr_s <- generated_data$gl_arr
    bed_arr_s <- generated_data$bed_arr

    true_surface_elevs <- generated_data$true_surface_elevs # May not use this... this is mostly for validation in state inference but we might not do that here
    true_thicknesses <- generated_data$true_thicknesses
    true_velocities <- generated_data$true_velocities

    # qsave(sim_results, file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))
    qsave(surface_obs_arr_s, file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
    qsave(friction_arr_s, file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
    qsave(gl_arr_s, file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))
    qsave(bed_arr_s, file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
    qsave(true_surface_elevs, file = paste0(train_data_dir, "/true_surface_elevs_", setf, "_", data_date, ".qs"))
    qsave(true_thicknesses, file = paste0(train_data_dir, "/true_thicknesses_", setf, "_", data_date, ".qs"))
    qsave(true_velocities, file = paste0(train_data_dir, "/true_velocities_", setf, "_", data_date, ".qs"))
}

t2 <- proc.time()
print(t2 - t1)

########################################
##          Plot simulations          ##
########################################

set <- sets[2]
setf <- formatC(set, width = 2, flag = "0")

### True thickness, friction, bed, grounding line
surface_obs_arr <- qread(file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))

friction_arr <- qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
bed_arr <- qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))

## Plot surface observations for a simulation
# png("plots/cnn/input/surface_obs_sim.png", width = 6, height = 8, units = "in", res = 300)
# par(mfrow = c(2, 1))
# matplot(surface_obs_arr[1, , , 1], type = "l", main = "Surface elevation (m)", ylab = "Surface elev. (m)", xlab = "Domain (grid points)")
# matplot(surface_obs_arr[1, , , 2], type = "l", main = "Velocity (m/yr)", ylab = "Velocity (m/yr)", xlab = "Domain (grid points)")
# dev.off()

## Check if the generated data is close to real data here
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs")
# vel_mat <- qread(file = "./data/velocity/all_velocity_arr.qs")
vel_mat <- qread(file = "./data/velocity/vel_smoothed.qs")

sim <- 1
png("./plots/cnn/input/sim_vs_real.png", width = 500, height = 600)
par(mfrow = c(2, 1))

matplot(surface_obs_arr[sim, , , 1],
    type = "l", col = "grey",
    xlab = "Grid point", ylab = "Surface elevation (m)"
)
matlines(surf_elev_mat, col = "salmon", lty = 2)
lines(ssa_steady$current_top_surface, col = "blue")
legend("topright", legend = c("Simulated", "Observed", "Steady state"), col = c("grey", "salmon", "blue"), lty = c(1, 2, 1), bty = "n")

matplot(vel_mat,
    ylim = c(0, 6000), type = "l", col = "salmon",
    xlab = "Grid point", ylab = "Velocity (m/a)"
)
# matlines(vel_mat2, col = "black", lty = 2)
matlines(surface_obs_arr[sim, , , 2], col = "grey", lty = 2)
lines(ssa_steady$current_velocity, col = "blue")
legend("topleft", legend = c("Observed", "Simulated", "Steady state"), col = c("salmon", "grey", "blue"), lty = c(1, 2, 1), bty = "n")
dev.off()

## Hovmoller plots of some simulations
plots <- list()

nsamples <- 4
sims <- sample(1:dim(surface_obs_arr)[1], size = nsamples)

space <- domain / 1000
time <- 1:years # 1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples * 4), nsamples, 4, byrow = T)

gl <- ceiling(gl_arr[1, 1] / (domain[length(domain)] / 1000) * length(domain))

# s <- 1
for (s in 1:nsamples) {
    sim <- sims[[s]]
    ind <- inds[s, ]
    surface_elev <- surface_obs_arr[sim, , , 1]
    velocity <- surface_obs_arr[sim, , , 2]
    grid_test$surface_elev <- as.vector(surface_elev)
    grid_test$velocity <- as.vector(velocity)

    surface_elev_plot <- ggplot(grid_test) +
        geom_tile(aes(space, time, fill = surface_elev)) +
        scale_y_reverse() +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        theme_bw() +
        theme(text = element_text(size = 24)) +
        # labs(fill="Thickness (m)")
        labs(fill = "Surface elev. (m)")

    velocity_plot <- ggplot(grid_test) +
        geom_tile(aes(space, time, fill = velocity)) +
        scale_y_reverse() +
        theme_bw() +
        theme(text = element_text(size = 24)) +
        scale_fill_distiller(palette = "Reds", direction = 1) +
        labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")"))


    # if (log_transform) {
    #   # fitted_fric_sim <- exp(fitted_friction[sim, 1:gl])
    #   friction_sim <- exp(friction_arr[sim, 1:gl])
    # } else {
    # fitted_fric_sim <- fitted_friction[sim, 1:gl]
    friction_sim <- friction_arr[sim, 1:gl]
    # }

    df <- data.frame(
        domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim # ,
        # fitted_fric = fitted_fric_sim
    )
    friction_plot <- ggplot(df, aes(x = domain, y = friction)) +
        geom_line() +
        # geom_line(aes(x = domain, y = fitted_fric), col = "red") +
        theme_bw() +
        xlim(c(0, domain[gl] / 1e3)) +
        xlab("Domain (km)") +
        ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")")) +
        theme(text = element_text(size = 24))

    bed_sim <- bed_arr[sim, ] #+ bed_mean
    #   fitted_bed_sim <- fitted_bed[sim, ] + bed_mean
    bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim) # , fitted_bed = fitted_bed_sim)

    bed_plot <- ggplot(bed_df, aes(x = domain, y = bed)) +
        geom_line() +
        # geom_line(aes(x = domain, y = fitted_bed), col = "red") +
        theme_bw() +
        ylim(c(-1500, -500)) +
        xlim(c(0, domain[gl] / 1e3)) +
        xlab("Domain (km)") +
        ylab("Bed (m)") +
        theme(text = element_text(size = 24))

    # plots[[ind[1]]] <- thickness_plot
    plots[[ind[1]]] <- surface_elev_plot
    plots[[ind[2]]] <- velocity_plot
    plots[[ind[3]]] <- friction_plot
    plots[[ind[4]]] <- bed_plot
}

png(
    file = paste0("./plots/cnn/input/simulations_", setf, "_", data_date, ".png"),
    width = 2400, height = 400 * nsamples
)
grid.arrange(grobs = plots, nrow = nsamples, ncol = 4)
dev.off()

# png(file = paste0("./plots/cnn/test.png"))
# matplot(surface_obs_arr[2,,,1], type = "l", col = "grey")
# matlines(surf_elev_mat, col = "red")
# dev.off()
