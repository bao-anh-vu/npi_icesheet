## Compute discrepancy between observed and simulated data
## for different friction coefficient fields

setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(dplyr)
library(fields)
library(Matrix)
library(matrixStats) # for the rowMaxs() function
library(R.utils)
library(mvtnorm)
library(abind)
library(parallel)
# library(mgcv)
library(FRK)

# source("./source/sim_params.R")
source("./source/sim_obs.R")
# source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/fit_basis.R")
source("./source/solve_ssa_nl_relax.R")
source("./source/surface_elev.R")
source("./source/solve_thickness.R")
source("./source/solve_velocity_azm.R")
source("./source/get_surface_obs.R")
# source("./source/azm_cond_sim.R")
source("./source/process_sim_results.R")

data_dir <- "./data/"
data_date <- "20241111" # "20241103"

## Flags
use_basal_melt_data <- T
leave_one_out <- T
nsims <- 100#0
avg_over_time <- F

# Physical params
# params <- qread(file = paste0("./data/training_data/phys_params_", data_date, ".qs"))
# params <- list(
#     secpera = 31556926, # seconds per annum
#     n = 3.0, # exponent in Glen's flow law
#     rho_i = 917.0, # ice density
#     rho_w = 1028.0, # sea water density
#     g = 9.81 # gravity constant
#     # A = 4.227e-25, #1.4579e-25, # flow rate parameter
# )

# params$m <- 1 / params$n
# params$B <- 0.7 * 1e6 * params$secpera^params$m
# params$A <- params$B^(-params$n)

params <- qread(file = paste0("./data/training_data/", "/phys_params_", data_date, ".qs"))


ssa_steady <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
domain <- ssa_steady$domain
J <- length(domain)

## Bed observations
bed_obs_df <- qread(file = paste0(data_dir, "bedmap/bed_obs_df_all.qs"))
bed_obs_chosen <- bed_obs_df
# bed_obs_df <- qread(file = paste0(data_dir, "/bed_obs_df.qs"))
# bed_obs_chosen <- bed_obs_df[bed_obs_df$chosen == 1, ]

## Read surface data
# vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_mat <- qread("./data/velocity/vel_smoothed.qs") # adjusted observed velocities
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

## Grounding line data
gl_pos <- qread(file = paste0(data_dir, "/grounding_line/gl_pos.qs"))
gl <- domain[gl_pos$ind] / 1e3
gl_ind <- gl_pos$ind

## SMB data
smb_data_racmo <- qread(file = paste0(data_dir, "/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
smb_avg <- colMeans(smb_data_racmo, na.rm = T)
params$as <- smb_avg # surface accumulation rate (m/s)

## Basal melt data (only available on shelves???)
if (use_basal_melt_data) {
    melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
    # qsave(flowline_shelf_melt, file = paste0(data_dir, "/SMB/flowline_shelf_melt.qs"))
    avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
    melt_nonmissing <- which(!is.na(avg_melt_rate))
    avg_melt_rate[1:(melt_nonmissing[1] - 1)] <- -0.1 # impose a melt rate of 1 m/a upstream of the first non-missing value
    avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) # for the remaining part of the shelf just use the last non-missing value
    avg_melt_rate <- -avg_melt_rate # inverting this as eventually smb is calculated as smb - melt
} else {
    avg_melt_rate <- c(rep(0.1, gl_ind), rep(3, J - gl_ind)) # simple basal melt model: 0.1 m/a acc on grounded ice, -0.5 m/a melt on shelf
}
params$ab <- avg_melt_rate # melt rate (m/s)

## Plot average melt rate
png(file = paste0("./plots/steady_state/avg_melt_rate_", data_date, ".png"), width = 800, height = 600)
plot(domain / 1000, params$ab,
    type = "l",
    xlab = "Distance along flowline (km)", ylab = "Basal melt rate (m/a)"
)
abline(v = ssa_steady$grounding_line[length(ssa_steady$grounding_line)], lty = 2, col = "red")
dev.off()

## Flowline data
# flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
# flowline <- qread(paste0(data_dir, "/flowline_regrid.qs"))
# J <- nrow(flowline) # number of grid points
# # flowline <- flowline[1:J, ]
# domain <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
# domain <- c(0, cumsum(na.omit(domain)))

######################################
##     Simulate bed and friction    ##
######################################

# sim_param_output <- sim_params(
#     nsims = nsims, domain = ssa_steady$domain,
#     bed_obs = bed_obs_df[bed_obs_df$chosen == 1, ]
# )

print("Simulating friction coefficient...")

fric_sims <- simulate_friction2(
    nsim = nsims, domain = domain
)

## Simulate beds
print("Simulating beds...")

# bed.sill <- 10e3
# bed.range <- 50e3
# bed.nugget <- 0 #200
# bed_sim_output <- simulate_bed(nsims, domain = domain,
#                         obs_locations = bed_obs_chosen$ind,
#                         obs = bed_obs_chosen$bed_elev) #,

# bed_sims <- bed_sim_output$sims
bed_prior <- qread(file = paste0("./data/bedmap/GP_fit_exp.qs"))
L <- t(chol(bed_prior$cov))
u_mat <- matrix(rnorm(nrow(L) * nsims), nrow = nrow(L), ncol = nsims)
bed_mean <- bed_prior$mean
mean_mat <- matrix(rep(bed_mean, nsims), nrow = nrow(bed_mean), ncol = nsims)
bed_sims <- mean_mat + L %*% u_mat # rnorm(nrow(L) * N)

## Need to fit basis here too
n_fric_basis <- 120
n_bed_basis <- 150

friction_basis <- fit_friction_basis(
    nbasis = n_fric_basis,
    domain = domain,
    fric_arr = t(fric_sims),
    log_transform = T,
    lengthscale = 3e3
)


bed_arr_demean <- bed_sims - mean_mat

## Fit basis to de-meaned bedrock
bed_basis <- fit_bed_basis(
    nbasis = n_bed_basis, domain = domain,
    bed_arr = t(bed_arr_demean),
    lengthscale = 2.5e3
)
bed_basis$mean <- bed_mean
bed_sims <- t(bed_basis$fitted_values) + mean_mat
# param_list <- lapply(1:nsims, function(r) {
#     list(
#         friction = fric_sims[, r], #* 1e6 * params$secpera^(1 / params$n),
#         bedrock = bed_sims[, r] #+ bed_mean
#     )

param_list <- lapply(1:nsims, function(r) {
    list(
        friction = exp(friction_basis$fitted_values[r, ]),
        bedrock = bed_basis$fitted_values[r, ] + bed_mean
    )
})

# sim_param_list <- sim_param_output$sim_param_list
# bed_sims <- sim_param_list$bedrock

# fric_scale <- 1e6 * params$secpera^(1 / params$n)
# fric_sims <- sim_param_list$friction #* fric_scale

## Run simulations for different friction coefficient fields
## starting from steady state
years <- dim(surf_elev_mat)[2] # number of years data is collected
warmup <- 0 # number of years to discard after steady state (for the model to adjust to new friction & bed)

sim_surface_obs <- list()
# for (s in 1:nsims) {
sim_out <- sim_obs(
    param_list = param_list,
    domain = domain,
    phys_params = params,
    years = years,
    warmup = warmup,
    # use_relaxation = T,
    ini_thickness = ssa_steady$current_thickness,
    ini_velocity = ssa_steady$current_velocity,
    smb = smb_avg,
    basal_melt = avg_melt_rate
    # log_transform = log_transform
)

prior_bad_sims <- sim_out$bad_sims
if (length(prior_bad_sims) > 0) {
    cat(length(prior_bad_sims), " prior simulations failed \n")
    sim_out$results <- sim_out$results[-prior_bad_sims]
    nsims <- length(sim_out$results)
}

generated_data <- process_sim_results(sims = sim_out$results)

sim_surface_obs <- generated_data$surface_obs_arr
# }

######################################
##     Compute model discrepancy    ##
######################################

if (leave_one_out) {
    years <- years - 1
    # vel_mat <- vel_mat[, 1:years]
    # surf_elev_mat <- surf_elev_mat[, 1:years]
}
se_sims <- lapply(1:nsims, function(s) sim_surface_obs[s, , 1:years, 1])
# se_sims <- lapply(sim_out, function(x) x$all_top_surface)
se_discr <- lapply(se_sims, function(M) surf_elev_mat[, 1:years] - M)

vel_sims <- lapply(1:nsims, function(s) sim_surface_obs[s, , 1:years, 2])
# vel_sims <- apply(sim_out, function(x) x$all_velocities)
vel_discr <- lapply(vel_sims, function(M) vel_mat[, 1:years] - M)

# GL_pos <- sapply(sim_out, function(x) x$grounding_line[length(x$grounding_line)])
GL_pos <- generated_data$gl_arr[, years]

pdf(file = paste0("./plots/discr/discrepancies_", data_date, ".pdf"), width = 10, height = 15)
nsims_plot <- 10
par(mfrow = c(2, 1))
matplot(domain / 1e3, bed_sims[, 1:nsims_plot],
    type = "l", col = "grey60", lwd = 2,
    xlab = "Distance along flowline",
    ylab = "Bed elevation (m)",
    main = ("Bed elevation simulations")
)
# lines(domain/1e3, ssa_steady$bedrock, col = "red", lwd = 2)
points(bed_obs_chosen$loc / 1e3, bed_obs_chosen$bed_elev, pch = 16, col = "salmon")

# png(file = paste0("./plots/discr/fric_sims_", data_date, ".png"), width = 800, height = 600)
matplot(domain / 1e3, fric_sims[, 1:nsims_plot],
    col = "grey60", type = "l", lwd = 2,
    xlab = "Distance along flowline",
    ylab = expression(paste("Friction coefficient (", m * a^{
        -1
    } * Pa^{
        -1 / 3
    }, ")")),
    main = ("Friction coefficient simulations")
)
# lines(domain/1e3, ssa_steady$friction_coef / fric_scale, col = "red", lwd = 2)

## Plot friction for different simulations

layout(matrix(1:nsims_plot, nrow = nsims_plot / 2, ncol = 2, byrow = TRUE))
for (s in 1:nsims_plot) {
    plot(domain / 1000, fric_sims[, s],
        type = "l", lty = 1, col = "salmon",
        # cex = 5,
        xlab = "Distance along flowline (km)",
        ylab = expression(paste("Friction coefficient (", m * a^{
            -1
        } * Pa^{
            -1 / 3
        }, ")")),
        main = paste0("Friction for simulation ", s)
    )
}

# Plot bed for different simulations
for (s in 1:nsims_plot) {
    plot(domain / 1000, bed_sims[, s],
        type = "l", lty = 1, col = "salmon",
        xlab = "Distance along flowline (km)", ylab = "Bed elevation (m)",
        main = paste0("Bed for simulation ", s)
    )
}

## Plot velocities for different simulations
# png(paste0("plots/discr/vel_sims_", data_date, ".png"), width = 1500, height = 2000)
# par(mfrow = c(nsims/2, 2))
# layout(matrix(1:nsims, nrow = nsims/2, ncol = 2, byrow = TRUE))

# Example: gradient from green to blue
n_lines <- ncol(vel_sims[[s]])  # number of lines to plot

# Create a color palette (green → blue)
cols <- cols <- adjustcolor(colorRampPalette(c("blue", "turquoise"))(n_lines), alpha.f = 0.6)

for (s in 1:nsims_plot) {
    matplot(domain / 1000, vel_mat,
        type = "l", lty = 1, col = "salmon",
        # cex = 5,
        ylim = c(0, 4000),
        xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
        main = paste0("Velocity for simulation ", s)
    )
    matlines(domain / 1000, vel_sims[[s]],
        type = "l", lty = 1, col = cols, lwd = 2
    )
    lines(domain/1000, ssa_steady$current_velocity, col = "black", lwd = 2)
    abline(v = GL_pos[s], lty = 2, col = "red")
    legend("bottomright", legend = c("Simulated", "Observed"), col = c("blue", "salmon"), lty = 1, bty = "n")
}

## Plot velocity discrepancy
# png(paste0("plots/discr/vel_discrepancy_", data_date, ".png"), width = 1500, height = 2000)
par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain / 1000, vel_discr[[s]],
        type = "l", lty = 1, # col = rgb(0,0,0,0.25),
        # cex = 5,
        col = cols, 
        ylim = c(-2000, 6000),
        xlab = "Distance along flowline (km)", ylab = "Velocity discrepancy (m/yr)",
        main = paste0("Velocity discrepancy for simulation ", s)
    )
    abline(h = 0, col = "grey", lty = 2)
    abline(v = GL_pos[s], lty = 2, col = "red")
}

## Plot surface elevation for different simulations
# png(paste0("plots/discr/se_sims_", data_date, ".png"), width = 1500, height = 2000)
par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain / 1000, surf_elev_mat,
        type = "l", lty = 1, col = "salmon",
        # cex = 5,
        ylim = c(0, 1600),
        xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
        main = paste0("Surface elevation for simulation ", s)
    )
    matlines(domain / 1000, se_sims[[s]],
        type = "l", lty = 1, col = cols, lwd = 2
    )
    lines(domain/1000, ssa_steady$current_top_surface, col = "black", lwd = 2)
    legend("topright", legend = c("Simulated", "Observed"), col = c("blue", "salmon"), lty = 1, bty = "n")
    abline(v = GL_pos[s], lty = 2, col = "red")
}

## Plot surface elevation discrepancy
# png(paste0("plots/discr/se_discrepancy_", data_date, ".png"), width = 1500, height = 2000)
par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain / 1000, se_discr[[s]],
        type = "l", lty = 1, # col = rgb(0,0,0,0.25),
        # cex = 5,
        col = cols,
        ylim = c(-500, 500),
        xlab = "Distance along flowline (km)", ylab = "Surface elevation discrepancy (m)",
        main = paste0("Surface elevation discrepancy for simulation ", s)
    )
    abline(h = 0, col = "grey", lty = 2)
    abline(v = GL_pos[s], lty = 2, col = "red")
}
dev.off()

years <- dim(vel_mat)[2]
    
if (avg_over_time) {
    ## Find average discrepancy (over simulations and over time)
    vel_discr_concat <- do.call(cbind, vel_discr)
    avg_vel_discr <- rowMeans(vel_discr_concat, na.rm = T)
    se_discr_concat <- do.call(cbind, se_discr)
    avg_se_discr <- rowMeans(se_discr_concat, na.rm = T)

    vel_discr_mat <- matrix(rep(avg_vel_discr, years), nrow = J, ncol = years)
    se_discr_mat <- matrix(rep(avg_se_discr, years), nrow = J, ncol = years)

} else {
    ## Find average discrepancy (over simulations) for each year
    vel_discr_mat <- Reduce("+", vel_discr) / length(vel_discr)
    se_discr_mat <- Reduce("+", se_discr) / length(se_discr)

    ## For the 11th year, just use the average discrepancy from the previous years
    avg_vel_discr <- rowMeans(vel_discr_mat, na.rm = T)
    avg_se_discr <- rowMeans(se_discr_mat, na.rm = T)
    vel_discr_mat <- cbind(vel_discr_mat, avg_vel_discr)
    se_discr_mat <- cbind(se_discr_mat, avg_se_discr)
}

adj_se_mat <- surf_elev_mat - se_discr_mat
adj_vel_mat <- vel_mat - vel_discr_mat

## Save discrepancy
if (avg_over_time) {
    file_tag <- "avg_"
} else {
    file_tag <- ""
}

qsave(vel_discr_mat, file = paste0(data_dir, "discrepancy/", file_tag, "vel_discr_", data_date, ".qs"))
qsave(se_discr_mat, file = paste0(data_dir, "discrepancy/", file_tag, "se_discr_", data_date, ".qs"))

# Save adjusted observed data
qsave(adj_se_mat, file = paste0(data_dir, "surface_elev/adj_se_mat_", file_tag, data_date, ".qs"))
qsave(adj_vel_mat, file = paste0(data_dir, "velocity/adj_vel_mat_", file_tag, data_date, ".qs"))

## Plot the average discrepancy
pdf(file = paste0("./plots/discr/adjusted_obs_", file_tag, data_date, ".pdf"), width = 15, height = 20, pointsize = 18)

par(mfrow = c(2, 1), cex = 1.5)
matplot(domain / 1000, vel_discr_mat,
    type = "l", col = cols, lwd = 2,
    # cex = 2,
    xlab = "Distance along flowline (km)", ylab = "Average velocity discrepancy (m/yr)",
    main = "Average velocity discrepancy"
)
abline(h = 0, col = "grey", lty = 2)
abline(v = gl, lty = 2)

matplot(domain / 1000, se_discr_mat,
    type = "l", col = cols, lwd = 2,
    # cex = 2,
    xlab = "Distance along flowline (km)", ylab = "Average surface elevation discrepancy (m)",
    main = "Average surface elevation discrepancy"
)
abline(h = 0, col = "grey", lty = 2)
abline(v = gl, lty = 2)

## Plot observed data minus discrepancy

# png(file = paste0("./plots/discr/vel_obs_minus_discr_", data_date, ".png"), width = 800, height = 600)
# pdf(file = paste0("./plots/discr/adjusted_obs_", data_date, ".pdf"), width = 10, height = 20)

par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain / 1000, vel_sims[[s]],
        type = "l", col = "black", lwd = 2,
        ylim = c(0, 4000),
        xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
        main = "(Adjusted) observed vs simulated velocity"
    )
    matlines(domain / 1000, vel_mat - vel_discr_mat,
        type = "l", col = "salmon", lwd = 2
    )
    legend("topright",
        legend = c("Simulated", "Observed - avg discrepancy"),
        col = c("grey", "salmon"), lwd = 2, bty = "n"
    )
    abline(v = gl, lty = 2)
}
# dev.off()

## Same plots for the surface elevation
# pdf(file = paste0("./plots/discr/adjusted_se_obs_", data_date, ".pdf"), width = 15, height = 20)


par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain / 1000, se_sims[[s]],
        type = "l", col = "black", lwd = 2,
        ylim = c(0, 1600),
        xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
        main = "(Adjusted) observed vs simulated surface elevation"
    )
    matlines(domain / 1000, adj_se_mat,
        type = "l", col = "salmon", lwd = 2
    )
    legend("topright",
        legend = c("Simulated", "Observed - avg discrepancy"),
        col = c("grey", "salmon"), lwd = 2, bty = "n"
    )
    abline(v = gl, lty = 2)
}
dev.off()
