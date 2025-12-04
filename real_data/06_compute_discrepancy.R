## Compute discrepancy between observed and simulated data
## for different friction coefficient fields

setwd("~/SSA_model/CNN/real_data/")

library(qs)
library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(tidyr)
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
# library(patchwork)

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
sets <- 51:100 #51:100 #51:100 # 51:100 #51:100 #6:20
setsf <- paste0("sets", sets[1], "-", sets[length(sets)])

## Flags
resimulate <- F
use_basal_melt_data <- T
leave_one_out <- T
nsims <- 100 # 0
use_relaxation <- F
warmup <- 5 # years to let the model "adjust" before actually collecting observations
avg_over_time <- T
fit_spline <- T

## Physical domain info
# ssa_steady <- qread(file = paste0(data_dir, "training_data/steady_state/steady_state_", data_date, ".qs"))
# domain <- ssa_steady$domain
# J <- length(domain)

flowline <- qread(paste0("./data/flowline_regrid.qs"))
J <- nrow(flowline) # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
domain <- c(0, cumsum(na.omit(flowline_dist)))

## Grounding line data
gl_pos <- qread(file = paste0(data_dir, "/grounding_line/gl_pos.qs"))
gl <- domain[gl_pos$ind] / 1e3
gl_ind <- gl_pos$ind

# Physical params

# params <- list(
#     secpera = 31556926, # seconds per annum
#     n = 3.0, # exponent in Glen's flow law
#     rho_i = 917.0, # ice density
#     rho_w = 1028.0, # sea water density
#     g = 9.81 # gravity constant
#     # A = 4.227e-25, #1.4579e-25, # flow rate parameter
# )

# params$m <- 1 / params$n
# params$B <- 0.6 * 1e6 * params$secpera^params$m
# params$A <- params$B^(-params$n)

# ## SMB data
# smb_data_racmo <- qread(file = paste0("./data/SMB/flowline_landice_smb.qs")) ## from 1979 to 2016
# smb_avg <- colMeans(smb_data_racmo, na.rm = T)
# params$as <- smb_avg # surface accumulation rate (m/s)

# if (use_basal_melt_data) {
#     melt_thwaites <- qread(file = "./data/SMB/flowline_shelf_melt.qs")
#     avg_melt_rate <- colMeans(melt_thwaites, na.rm = T)
#     # melt_nonmissing <- which(!is.na(avg_melt_rate))
#     # avg_melt_rate[1:(melt_nonmissing[1] - 1)] <- -0.1 # seq(0, avg_melt_rate[melt_nonmissing[1]], length.out = melt_nonmissing[1]-1)
#     # avg_melt_rate[is.na(avg_melt_rate)] <- tail(avg_melt_rate[melt_nonmissing], 1) # mean(avg_melt_rate[melt_nonmissing])
#     # avg_melt_rate <- -avg_melt_rate # inverting this as eventually smb is calculated as smb - melt

#     avg_melt_shelf <- mean(melt_thwaites, na.rm = T)
#     # avg_melt_rate <- rep(-0.1, J)
#     avg_melt_rate <- rep(0, J)
#     avg_melt_rate[gl_ind:J] <- -avg_melt_shelf # extra minus sign here as the melt rate in my model has the opposite sign as what is in the data
# } else { ## assume no melt
#     avg_melt_rate <- rep(0, J)
# }
# params$ab <- avg_melt_rate # melt rate (m/s)

params <- qread(file = paste0("./data/training_data/phys_params_", data_date, ".qs"))

## Bed observations
bed_obs_df <- qread(file = paste0(data_dir, "bedmap/bed_obs_df_all.qs"))
bed_obs_chosen <- bed_obs_df

## Read surface data
# vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_mat <- qread("./data/velocity/vel_smoothed.qs") # adjusted observed velocities
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs") # this is on grounded ice only

## Initial velocity field
vel_curr_smooth <- qread(file = paste0("./data/velocity/ini_vel_", data_date, ".qs"))
vel_err_sd <- qread(file = paste0("./data/velocity/vel_err_sd_", data_date, ".qs"))

# ## Plot average melt rate
# png(file = paste0("./plots/steady_state/avg_melt_rate_", data_date, ".png"), width = 800, height = 600)
# plot(domain / 1000, params$ab,
#     type = "l",
#     xlab = "Distance along flowline (km)", ylab = "Basal melt rate (m/a)"
# )
# abline(v = ssa_steady$grounding_line[length(ssa_steady$grounding_line)], lty = 2, col = "red")
# dev.off()


######################################
##     Simulate bed and friction    ##
######################################

# sim_param_output <- sim_params(
#     nsims = nsims, domain = ssa_steady$domain,
#     bed_obs = bed_obs_df[bed_obs_df$chosen == 1, ]
# )

if (resimulate) {
    print("Simulating friction coefficient...")
    set.seed(2025)
    fric_sims <- simulate_friction2(
    nsim = nsims, domain = domain
)

## Simulate beds
print("Simulating beds...")
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

}



# sim_param_list <- sim_param_output$sim_param_list
# bed_sims <- sim_param_list$bedrock

# fric_scale <- 1e6 * params$secpera^(1 / params$n)
# fric_sims <- sim_param_list$friction #* fric_scale

years <- dim(surf_elev_mat)[2] # number of years data is collected
# warmup <- 0 # number of years to discard after steady state (for the model to adjust to new friction & bed)

######################################
##     Compute model discrepancy    ##
######################################
if (resimulate) {
    ## Run simulations for different friction coefficient fields
    ## starting from steady state

    sim_surface_obs <- list()
    # for (s in 1:nsims) {
    sim_out <- sim_obs(
        param_list = param_list,
        domain = domain,
        phys_params = params,
        years = years,
        warmup = warmup,
        use_relaxation = use_relaxation,
        relax_years = warmup, # over how many years to relax towards observed thickness
        # ini_thickness = ssa_steady$current_thickness,
        # ini_surface = ssa_steady$current_top_surface,
        # ini_velocity = ssa_steady$current_velocity,
        ini_surface = surf_elev_mat[, 1],
        ini_velocity = vel_curr_smooth,
        vel_err_sd = vel_err_sd
        # smb = smb_avg,
        # basal_melt = avg_melt_rate
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
    
    se_sims <- lapply(1:nsims, function(s) sim_surface_obs[s, , 1:years, 1])
    # se_sims <- lapply(sim_out, function(x) x$all_top_surface)

    vel_sims <- lapply(1:nsims, function(s) sim_surface_obs[s, , 1:years, 2])
    # vel_sims <- apply(sim_out, function(x) x$all_velocities)

    qsave(bed_sims, file = paste0(data_dir, "discrepancy/", setsf, "/bed_sims_", data_date, ".qs"))
    qsave(fric_sims, file = paste0(data_dir, "discrepancy/", setsf, "/fric_sims_", data_date, ".qs"))
    qsave(vel_sims, file = paste0(data_dir, "discrepancy/", setsf, "/vel_sims_list_", data_date, ".qs"))
    qsave(se_sims, file = paste0(data_dir, "discrepancy/", setsf, "/se_sims_list_", data_date, ".qs"))
} else {
    bed_sims <- qread(file = paste0(data_dir, "discrepancy/", setsf, "/bed_sims_", data_date, ".qs"))
    fric_sims <- qread(file = paste0(data_dir, "discrepancy/", setsf, "/fric_sims_", data_date, ".qs"))
    vel_sims <- qread(file = paste0(data_dir, "discrepancy/", setsf, "/vel_sims_list_", data_date, ".qs"))
    se_sims <- qread(file = paste0(data_dir, "discrepancy/", setsf, "/se_sims_list_", data_date, ".qs"))
}

if (leave_one_out) {
    years <- years - 1
    vel_mat <- vel_mat[, 1:years]
    surf_elev_mat <- surf_elev_mat[, 1:years]
}
se_discr <- lapply(se_sims, function(M) surf_elev_mat[, 1:years] - M[, 1:years])
vel_discr <- lapply(vel_sims, function(M) vel_mat[, 1:years] - M[, 1:years])


# GL_pos <- sapply(sim_out, function(x) x$grounding_line[length(x$grounding_line)])
# GL_pos <- generated_data$gl_arr[, years]

### Plot model discrepancy
nsims_plot <- 4
plot_range <- 1:sum(!is.na(surf_elev_mat[, 1]))

# pdf(file = paste0("./plots/discr/discrepancies_", data_date, ".pdf"), width = 10, height = 15)
png(file = paste0("./plots/discr/bed_fric_sims_", data_date, ".png"), width = 1500, height = 250 * nsims_plot, res = 200)

par(mfrow = c(2, 1))
matplot(domain / 1e3, bed_sims[, 1:nsims_plot],
    type = "l", col = "grey60", lwd = 1.5,
    xlab = "Distance along flowline",
    ylab = "Bed elevation (m)",
    main = ("Bed elevation simulations")
)
# lines(domain/1e3, ssa_steady$bedrock, col = "red", lwd = 1.5)
points(bed_obs_chosen$loc / 1e3, bed_obs_chosen$bed_elev, pch = 16, col = "salmon")

# png(file = paste0("./plots/discr/fric_sims_", data_date, ".png"), width = 800, height = 600)
matplot(domain / 1e3, fric_sims[, 1:nsims_plot],
    col = "grey60", type = "l", lwd = 1.5,
    xlab = "Distance along flowline",
    ylab = expression(paste("Friction coefficient (", m * a^{
        -1
    } * Pa^{
        -1 / 3
    }, ")")),
    main = ("Friction coefficient simulations")
)
# lines(domain/1e3, ssa_steady$friction_coef / fric_scale, col = "red", lwd = 1.5)
dev.off()

## Plot friction for different simulations
png(paste0("plots/discr/fric_sims_", data_date, ".png"), width = 1500, height = 250 * nsims_plot)
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
dev.off()

# Plot bed for different simulations
png(paste0("plots/discr/bed_sims_", data_date, ".png"), width = 1500, height = 300 * nsims_plot)
for (s in 1:nsims_plot) {
    plot(domain / 1000, bed_sims[, s],
        type = "l", lty = 1, col = "salmon",
        xlab = "Distance along flowline (km)", ylab = "Bed elevation (m)",
        main = paste0("Bed for simulation ", s)
    )
}
dev.off()

## Plot velocities for different simulations
png(paste0("plots/discr/vel_sims_", data_date, ".png"), width = 1500, height = 300 * nsims_plot, res = 200)
par(mfrow = c(nsims_plot / 2, 2))
# layout(matrix(1:nsims, nrow = nsims/2, ncol = 2, byrow = TRUE))

# Example: gradient from green to blue
n_lines <- years+1 # number of lines to plot
plot_range <- 1:gl_ind

# Create a color palette (green → blue)
cols <- adjustcolor(colorRampPalette(c("maroon", "lightpink"))(n_lines), alpha.f = 0.8)
# grey_cols <- adjustcolor(colorRampPalette(c("lightcoral", "mistyrose"))(n_lines), alpha.f = 0.8)
grey_cols <- adjustcolor(colorRampPalette(c("grey40", "grey80"))(n_lines), alpha.f = 0.8)

for (s in 1:nsims_plot) {
    matplot(domain / 1000, vel_mat,
        type = "l", lty = 1, col = "salmon",
        # cex = 5,
        ylim = c(0, 4000),
        xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
        main = paste0("Velocity for simulation ", s)
    )
    matlines(domain / 1000, vel_sims[[s]],
        type = "l", lty = 1, col = cols, lwd = 1.5
    )
    # lines(domain / 1000, ssa_steady$current_velocity, col = "black", lwd = 1.5)
    # abline(v = GL_pos[s], lty = 2, col = "red")
    legend("bottomright", legend = c("Simulated", "Observed"), col = c("blue", "salmon"), lty = 1, bty = "n")
}
dev.off()


## Plot velocity discrepancy
png(paste0("plots/discr/vel_discrepancy_", data_date, ".png"), width = 1500, height = 300 * nsims_plot, res = 200)
par(mfrow = c(nsims_plot / 2, 2))
for (s in 1:nsims_plot) {
    matplot(domain[plot_range] / 1000, vel_discr[[s]][plot_range, ],
        type = "l", lty = 1, # col = rgb(0,0,0,0.25),
        # cex = 5,
        col = cols,
        ylim = c(-1000, 1500),
        xlab = "Distance along flowline (km)", ylab = "Discrepancy (m/yr)",
        main = paste0("Simulation ", s)
    )
    abline(h = 0, col = "grey", lty = 2)
    # abline(v = GL_pos[s], lty = 2, col = "black")
}
dev.off()

## Plot surface elevation for different simulations
png(paste0("plots/discr/se_sims_", data_date, ".png"), width = 1500, height = 550 * nsims_plot, res = 300)
par(mfrow = c(nsims_plot, 2))
for (s in 1:nsims_plot) {
    matplot(domain[plot_range] / 1000, se_sims[[s]][plot_range, ],
        type = "l", lty = 1, lwd = 1.5, col = grey_cols,
        # cex = 5,
        ylim = c(0, 1400),
        xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
        main = paste0("Surface elevation for simulation ", s)
    )
    matlines(domain[plot_range] / 1000, surf_elev_mat[plot_range, ],
        type = "l", lty = 1, col = cols, lwd = 1.5
    )
    # lines(domain / 1000, ssa_steady$current_top_surface, col = "black", lwd = 1.5)
    # legend("topright", legend = c("Simulated", "Observed"), col = c("blue", "salmon"), lty = 1, bty = "n")
    # abline(v = GL_pos[s], lty = 2, col = "red")

    matplot(domain[plot_range] / 1000, vel_sims[[s]][plot_range, ],
        type = "l", lty = 1, lwd = 1.5, col = grey_cols,
        # cex = 5,
        ylim = c(0, 2600),
        xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
        main = paste0("Velocity for simulation ", s)
    )
    matlines(domain[plot_range] / 1000, vel_mat[plot_range, ],
        type = "l", lty = 1, col = cols, lwd = 1.5
    )
    # lines(domain / 1000, ssa_steady$current_velocity, col = "black", lwd = 1.5)
    # abline(v = GL_pos[s], lty = 2, col = "red")
    # legend("bottomright", legend = c("Simulated", "Observed"), col = c("blue", "salmon"), lty = 1, bty = "n")

}
dev.off()

## Plot surface elevation discrepancy
# png(paste0("plots/discr/se_discrepancy_", data_date, ".png"), width = 1500, height = 550 * nsims_plot, res = 300)
# par(mfrow = c(nsims_plot, 2))
# for (s in 1:nsims_plot) {
#     matplot(domain[plot_range] / 1000, se_discr[[s]][plot_range, ],
#         type = "l", lty = 1, # col = rgb(0,0,0,0.25),
#         cex = 1.5,
#         col = cols,
#         # ylim = c(-500, 500),
#         xlab = "Distance along flowline (km)", ylab = "Surface elev. discrepancy (m)",
#         main = paste0("Simulation ", s)
#     )
#     abline(h = 0, col = "grey", lty = 2)
#     # abline(v = GL_pos[s], lty = 2, col = "red")

#     matplot(domain[plot_range] / 1000, vel_discr[[s]][plot_range, ],
#         type = "l", lty = 1, # col = rgb(0,0,0,0.25),
#         cex = 1.5,
#         col = cols,
#         ylim = c(-1000, 1500),
#         xlab = "Distance along flowline (km)", ylab = "Velocity discrepancy (m/yr)",
#         main = paste0("Simulation ", s)
#     )
#     abline(h = 0, col = "grey", lty = 2)
# }
# dev.off()

## ggplot equivalent
x_vals <- domain[plot_range] / 1000

tidy_mat <- function(mat, sim_id, varname) {
  as.data.frame(mat[plot_range, ]) %>%
    mutate(x = x_vals) %>%
    pivot_longer(
      cols = -x,
      names_to = "line",
      values_to = "value"
    ) %>%
    mutate(sim = sim_id, variable = varname)
}

## Plot simulations 
se_sim_plots <- list()
vel_sim_plots <- list()

for (s in 1:nsims_plot) {

  # --- Build tidy data for sim s ---
  df_surf_sim <- tidy_mat(se_sims[[s]][, 1:years], s, "Surface elevation")
  df_surf_obs <- tidy_mat((surf_elev_mat)[, 1:years], s, "Surface elevation")

  df_vel_sim  <- tidy_mat(vel_sims[[s]][, 1:years], s, "Velocity")
  df_vel_obs  <- tidy_mat((vel_mat)[, 1:years], s, "Velocity")

  # --- Surface elevation panel ---
  p1 <- ggplot() +
    # simulated curves (grey)
    geom_line(
      data = df_surf_sim,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = grey_cols, guide = "none") +
    new_scale_color() +
    # observed minus discrepancy curves (colored)
    geom_line(
      data = df_surf_obs,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = cols, guide = "none") +

    coord_cartesian(ylim = c(0, 1400)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Surface elevation (m)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # --- Velocity panel ---
  p2 <- ggplot() +
    # simulated curves (grey)
    geom_line(
      data = df_vel_sim,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = grey_cols, guide = "none") +
    new_scale_color() +
    # observed minus discrepancy curves (colored)
    geom_line(
      data = df_vel_obs,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = cols, guide = "none") +

    coord_cartesian(ylim = c(0, 2800)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Velocity (m/yr)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # Save row of two plots
    se_sim_plots[[s]] <- p1 
    vel_sim_plots[[s]] <- p2
}

png(paste0("plots/discr/sims_ggplot_", data_date, ".png"), width = 1800, height = 600 * nsims_plot, res = 300)
grid.arrange(grobs = c(se_sim_plots, vel_sim_plots), 
            layout_matrix = matrix(1:(nsims_plot * 2), nsims_plot, 2) )
dev.off()


## Plot discrepancy using ggplot

se_discr_plots <- list()
vel_discr_plots <- list()

for (s in 1:nsims_plot) {

  # Tidy data for sim s
  df_se  <- tidy_mat(se_discr[[s]],  s, "Surface elev. discrepancy")
  df_vel <- tidy_mat(vel_discr[[s]], s, "Velocity discrepancy")

  # --- Surface Elevation Discrepancy Plot ---
  p1 <- ggplot(df_se, aes(x = x, y = value, group = line, color = line)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(legend.position = "none", 
        text = element_text(size = 14)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Surface elevation \n discrepancy (m)"
    ) 

  # --- Velocity Discrepancy Plot ---
  p2 <- ggplot(df_vel, aes(x = x, y = value, group = line, color = line)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = c(-1000, 1500)) +
    theme_bw() +
    theme(legend.position = "none", 
        text = element_text(size = 14)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Velocity \n discrepancy (m/yr)"
    ) 

  # store combined 2-panel row
  se_discr_plots[[s]] <- p1
  vel_discr_plots[[s]] <- p2
}

se_discr_plots <- lapply(se_discr_plots, function(p) p + theme(plot.margin = margin(20, 20, 20, 20)))
vel_discr_plots <- lapply(vel_discr_plots, function(p) p + theme(plot.margin = margin(20, 20, 20, 20)))

png(paste0("plots/discr/discrepancies_ggplot_", data_date, ".png"), width = 2000, height = 700 * nsims_plot, res = 200)
grid.arrange(grobs = c(se_discr_plots, vel_discr_plots), 
            layout_matrix = matrix(1:(nsims_plot * 2), nsims_plot, 2))
# print(plots[[1]])
dev.off()


## Find average discrepancy (over simulations) for each year
# years <- dim(vel_mat)[2]
vel_discr_mat <- Reduce("+", vel_discr) / length(vel_discr)
se_discr_mat <- Reduce("+", se_discr) / length(se_discr)

# vel_discr_sub <- vel_discr[1:100]
# se_discr_sub <- se_discr[1:100]

# vel_discr_mat_sub <- Reduce("+", vel_discr_sub) / length(vel_discr_sub)
# se_discr_mat_sub <- Reduce("+", se_discr_sub) / length(se_discr_sub)


if (avg_over_time) {
    avg_se_discr <- rowMeans(se_discr_mat, na.rm = T)
    avg_vel_discr <- rowMeans(vel_discr_mat, na.rm = T)

    # avg_se_discr_sub <- rowMeans(se_discr_mat_sub, na.rm = T)
    # avg_vel_discr_sub <- rowMeans(vel_discr_mat_sub, na.rm = T)

    if (fit_spline) {
        ## Fit a spline through the avg discrepancy
        library(mgcv)
        spline_fit <- gam(avg_se_discr ~ s(domain, k = 20))
        se_discr_smooth <- predict(spline_fit, newdata = data.frame(domain = domain))
        # se_discr_smooth <- se_discr_smooth[!is.na(avg_se_discr)]

        spline_fit <- gam(avg_vel_discr ~ s(domain, k = 20))
        vel_discr_smooth <- predict(spline_fit, newdata = data.frame(domain = domain))

        png(paste0("./plots/discr/se_discrepancy_smooth.png"), width = 1500, height = 500, res = 200)
        par(mfrow = c(1, 2))
        plot(domain[plot_range] / 1e3, avg_se_discr[plot_range],
            type = "l", col = cols[5], lwd = 1.5, lty = 1,
            xlab = "Distance along flowline (km)", ylab = "Discrepancy (m)",
            main = "Surface elevation discrepancy"
        )
        lines(domain[plot_range] / 1e3, se_discr_smooth[plot_range], col = "black", lwd = 1.5)

        plot(domain[plot_range] / 1e3, avg_vel_discr[plot_range],
            type = "l", col = cols[5], lwd = 1.5, lty = 1,
            xlab = "Distance along flowline (km)", ylab = "Discrepancy (m/yr)",
            main = "Velocity discrepancy"
        )
        lines(domain[plot_range] / 1e3, vel_discr_smooth[plot_range], col = "black", lwd = 1.5)
        dev.off()

        avg_se_discr[!is.na(avg_se_discr)] <- se_discr_smooth[!is.na(avg_se_discr)]
        avg_vel_discr[!is.na(avg_vel_discr)] <- vel_discr_smooth[!is.na(avg_vel_discr)]
        # se_discr_smooth[is.na(avg_se_discr)] <- NA
        # vel_discr_smooth[is.na(avg_vel_discr)] <- NA
        # qsave(se_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/se_discr_smooth_", data_date, ".qs"))
        # qsave(vel_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/vel_discr_smooth_", data_date, ".qs"))
    }

    se_discr_mat <- matrix(rep(avg_se_discr, years), nrow = length(se_discr_smooth), ncol = years)
    vel_discr_mat <- matrix(rep(avg_vel_discr, years), nrow = length(vel_discr_smooth), ncol = years)

    qsave(se_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/se_discr_avg_", data_date, ".qs"))
    qsave(vel_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/vel_discr_avg_", data_date, ".qs"))
} else { # save the non-time-averaged discrepancy matrices

    qsave(vel_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/vel_discr_", data_date, ".qs"))
    qsave(se_discr_mat, file = paste0(data_dir, "discrepancy/", setsf, "/se_discr_", data_date, ".qs"))
}

# # if (avg_over_time) {
#     # ## Find average discrepancy (over simulations and over time)
#     # vel_discr_concat <- do.call(cbind, vel_discr)
#     # avg_vel_discr <- rowMeans(vel_discr_concat, na.rm = T)
#     # se_discr_concat <- do.call(cbind, se_discr)
#     # avg_se_discr <- rowMeans(se_discr_concat, na.rm = T)

#     # vel_discr_mat <- matrix(rep(avg_vel_discr, years), nrow = J, ncol = years)
#     # se_discr_mat <- matrix(rep(avg_se_discr, years), nrow = J, ncol = years)

# # } else {
#     ## Find average discrepancy (over simulations) for each year
#     vel_discr_mat <- Reduce("+", vel_discr) / length(vel_discr)
#     se_discr_mat <- Reduce("+", se_discr) / length(se_discr)

#     # ## For the 11th year, just use the average discrepancy from the previous years
#     # avg_vel_discr <- rowMeans(vel_discr_mat, na.rm = T)
#     # avg_se_discr <- rowMeans(se_discr_mat, na.rm = T)
#     # vel_discr_mat <- cbind(vel_discr_mat, avg_vel_discr)
#     # se_discr_mat <- cbind(se_discr_mat, avg_se_discr)
# # }



#################################################

## Plot the average discrepancy
# pdf(file = paste0("./plots/discr/adjusted_obs_", data_date, ".pdf"), width = 15, height = 20, pointsize = 18)
png(file = paste0("./plots/discr/avg_discrepancy_", data_date, ".png"), width = 1500, height = 2000, res = 300)
par(mfrow = c(2, 1))
matplot(domain / 1000, se_discr_mat,
    type = "l", col = cols, lwd = 1.5, lty = 1,
    # cex = 2,
    xlab = "Distance along flowline (km)", ylab = "Discrepancy (m)",
    main = "Average surface elevation discrepancy"
)
abline(h = 0, col = "grey", lty = 2)
abline(v = gl, lty = 2)

matplot(domain / 1000, vel_discr_mat,
    type = "l", col = cols, lwd = 1.5, lty = 1,
    # cex = 2,
    xlab = "Distance along flowline (km)", ylab = "Discrepancy (m/yr)",
    main = "Average velocity discrepancy"
)
abline(h = 0, col = "grey", lty = 2)
abline(v = gl, lty = 2)
dev.off()

# Plot observed data minus discrepancy

# pdf(file = paste0("./plots/discr/adjusted_obs_", data_date, ".pdf"), width = 10, height = 20)
# png(file = paste0("./plots/discr/adjusted_obs_", data_date, ".png"), width = 1500, height = 550 * nsims_plot, res = 300)
# par(mfrow = c(nsims_plot, 2))
# for (s in 1:nsims_plot) {
#     # legend("topright",
#     #     legend = c("Simulated", "Observed - avg discrepancy"),
#     #     col = c("black", "salmon"), lwd = 1.5, bty = "n"
#     # )
#     # abline(v = gl, lty = 2)

#     matplot(domain[plot_range] / 1000, se_sims[[s]][plot_range, ],
#         type = "l", col = grey_cols, lwd = 1.5, lty = 1,
#         ylim = c(0, 1400),
#         xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)",
#         main = paste0("Simulation ", s)
#     )
#     matlines(domain[plot_range] / 1000, (surf_elev_mat - se_discr_mat)[plot_range, ],
#         type = "l", col = cols, lwd = 1.5, lty = 1
#     )
#     # legend("topright",
#     #     legend = c("Simulated", "Observed - avg discrepancy"),
#     #     col = c("black", "salmon"), lwd = 1.5, bty = "n"
#     # )
#     # abline(v = gl, lty = 2)

#     matplot(domain[plot_range] / 1000, vel_sims[[s]][plot_range, ],
#         type = "l", col = grey_cols, lwd = 1.5, lty = 1,
#         ylim = c(0, 3000),
#         xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)",
#         main = paste0("Simulation ", s)
#     )
#     matlines(domain[plot_range] / 1000, (vel_mat - vel_discr_mat)[plot_range, ],
#         type = "l", col = cols, lwd = 1.5, lty = 1
#     )
# }
# dev.off()

## Same plot but in ggplot

adj_se_sim_plots <- list()
adj_vel_sim_plots <- list()

for (s in 1:nsims_plot) {

  # --- Build tidy data for sim s ---
  df_surf_sim <- tidy_mat(se_sims[[s]][, 1:years], s, "Surface elevation")
  df_surf_obs <- tidy_mat((surf_elev_mat - se_discr_mat)[, 1:years], s, "Surface elevation")

  df_vel_sim  <- tidy_mat(vel_sims[[s]][, 1:years], s, "Velocity")
  df_vel_obs  <- tidy_mat((vel_mat - vel_discr_mat)[, 1:years], s, "Velocity")

  # --- Surface elevation panel ---
  p1 <- ggplot() +
    # simulated curves (grey)
    geom_line(
      data = df_surf_sim,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = grey_cols, guide = "none") +
    new_scale_color() +
    # observed minus discrepancy curves (colored)
    geom_line(
      data = df_surf_obs,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = cols, guide = "none") +

    coord_cartesian(ylim = c(0, 1400)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Surface elevation (m)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # --- Velocity panel ---
  p2 <- ggplot() +
    # simulated curves (grey)
    geom_line(
      data = df_vel_sim,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = grey_cols, guide = "none") +
    new_scale_color() +
    # observed minus discrepancy curves (colored)
    geom_line(
      data = df_vel_obs,
      aes(x = x, y = value, group = line, color = line),
      linewidth = 0.5
    ) +
    scale_color_manual(values = cols, guide = "none") +

    coord_cartesian(ylim = c(0, 2800)) +
    labs(
      title = paste("Simulation", s),
      x = "Distance along flowline (km)",
      y = "Velocity (m/yr)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  # Save row of two plots
    adj_se_sim_plots[[s]] <- p1 
    adj_vel_sim_plots[[s]] <- p2
}

png(paste0("plots/discr/adjusted_obs_ggplot_", data_date, ".png"), width = 2000, height = 700 * nsims_plot, res = 300)
grid.arrange(grobs = c(adj_se_sim_plots, adj_vel_sim_plots), 
            layout_matrix = matrix(1:(nsims_plot * 2), nsims_plot, 2) )
dev.off()

# ## Plot adjusted obs by year
# se_sims_yr <- lapply(se_sims, function(M) M[, yr])
# vel_sims_yr <- lapply(vel_sims, function(M) M[, yr])
# se_sims_mat <- do.call(cbind, se_sims_yr)
# vel_sims_mat <- do.call(cbind, vel_sims_yr)
# nsims_plot <- 10
# png(file = paste0("./plots/discr/adjusted_obs_by_year_", data_date, ".png"), 
#         width = 1500, height = 550 * nsims_plot, res = 300)
#         par(mfrow = c(years, 2))

# for (yr in 1:years) {
    
#         matplot(domain[plot_range] / 1000, se_sims_mat[plot_range, 1:nsims_plot],
#             type = "l", col = "grey", lwd = 1.5, lty = 1,
#             ylim = c(0, 1600),
#             xlab = "Distance along flowline (km)", ylab = "Surface elevation (m)") #,
#             # main = paste0("(Adjusted) observed vs simulated surface elevation - Year ", yr)
        
#         matlines(domain[plot_range] / 1000, (surf_elev_mat - se_discr_mat)[plot_range, yr],
#             type = "l", col = "salmon", lwd = 1.5, lty = 1
#         )

#         matplot(domain[plot_range] / 1000, vel_sims_mat[plot_range, 1:nsims_plot],
#             type = "l", col = "grey", lwd = 1.5, lty = 1,
#             ylim = c(0, 4000),
#             xlab = "Distance along flowline (km)", ylab = "Velocity (m/yr)") #,
#             # main = paste0("(Adjusted) observed vs simulated velocity - Year ", yr)
        
#         matlines(domain[plot_range] / 1000, (vel_mat - vel_discr_mat)[plot_range, yr],
#             type = "l", col = "salmon", lwd = 1.5, lty = 1
#         )
    
# }
# dev.off()
#_____________________
## From each simulation, extract the first year in se_sims and vel_sims
# se_sims_ls <- list()
# vel_sims_ls <- list()
# for (yr in 1:years) {

#     ## Turn into dataframe with year column
#     se_sims_df <- do.call(cbind, se_sims_yr)
#     vel_sims_df <- do.call(cbind, vel_sims_yr)
#     se_sims_df <- as.data.frame(se_sims_df)
#     vel_sims_df <- as.data.frame(vel_sims_df)
#     se_sims_df %>% mutate(
#         distance = domain,
#         year = yr
#     ) 
#     vel_sims_df %>% mutate(
#         distance = domain,
#         year = yr
#     )
#     se_sims_ls[[yr]] <- se_sims_df
#     vel_sims_ls[[yr]] <- vel_sims_df
# }

# se_sims_df <- bind_rows(se_sims_ls)
# vel_sims_df <- bind_rows(vel_sims_ls)

# ## Do panel plot of adjusted observed vs simulated for each year
# # library(ggplot2)
# # library(tidyr)
# library(reshape2)
# se_sims_melt <- melt(se_sims_df, id.vars = c("distance", "year"), variable.name = "simulation", value.name = "surface_elev")
# vel_sims_melt <- melt(vel_sims_df, id.vars = c("distance", "year"), variable.name = "simulation", value.name = "velocity")

## Plot simulated vs observed - discrepancy in one plot
# png(file = paste0("./plots/discr/first_year_adjusted_obs_", data_date, ".png"), 
#     width = 1500, height = 550 * nsims_plot, res = 300)

