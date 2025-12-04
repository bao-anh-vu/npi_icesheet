## Generate 1000 simulations

setwd("~/SSA_model/CNN/simbed")

rm(list = ls())

library(parallel)
library(Matrix)
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
library(R.utils)
# library("sp")
library(fields)
# library("tidyr")
library(dplyr)
library(matrixStats) # for the rowMaxs() function
library(mvtnorm)
library(abind)
library(ggplot2)
# library("plotly")
library(gridExtra)
library(FRK)
library(qs)

source("./source/sim_params.R")
# source("./source/run_sims.R")
source("./source/sim_obs.R")
source("./source/process_sim_results.R")
source("./source/fit_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
source("./source/create_ref.R")
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/get_surface_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")

## Some flags
regenerate_sims <- F
reprocess_sim_results <- F
# refit_basis <- T
save_sims <- F
log_transform <- T
# sim_beds <- T

train_data_dir <- "./training_data"

## Presets
data_date <- "20240320" # "20220329"
N <- 1000 # number of simulations per set
warmup <- 0
years <- 20 + warmup
sets <- 1#:10 #11:20
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 150

# 0. Load ice sheet at steady state
ssa_steady <- readRDS(file = paste0(train_data_dir, "/initial_conds/ssa_steady_20220329.rds", sep = ""))
domain <- ssa_steady$domain

# png(paste0("./plots/temp/steady_state.png"), width = 800, height = 800)
# plot(domain, ssa_steady$top_surface, type = 'l', ylim = c(-2000, 2000))
# lines(domain, ssa_steady$bottom_surface)
# lines(domain, ssa_steady$bedrock)
# dev.off()

# 0. Simulate bed observations
set.seed(2024)
ref_bed <- ssa_steady$bedrock
## Randomly select 50 locations between x = 0 and x = 800 km
n_bed_obs <- 50
obs_ind <- sort(sample(length(ref_bed), n_bed_obs))
obs_bed <- ref_bed[obs_ind] + rnorm(n_bed_obs, mean = 0, sd = 20) ## add noise
bed_obs <- list(locations = obs_ind, obs = obs_bed)
rm(.Random.seed, envir = globalenv())

if (save_sims) {
  qsave(bed_obs, file = paste0(train_data_dir, "/bed_obs_", data_date, ".qs"))
}

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1 / 3)

sim_results_list <- list()
flags <- c()

if (regenerate_sims) {

  for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Generating set", set, "\n")

    sim_param_list <- sim_params(
      nsims = N, domain = ssa_steady$domain,
      bed_obs = bed_obs
    )

    ## Then fit basis here
    cat("Fitting basis functions for set", set, "\n")
    setf <- formatC(set, width = 2, flag = "0")

    ## Fit basis to log(friction)
    friction_basis <- fit_friction_basis(
      nbasis = 150,
      domain = domain,
      fric_arr = sim_param_list$friction,
      log_transform = log_transform, 
      lengthscale = 8e3
    )



    png(file = paste0("./plots/temp/friction_basis_", setf, "_", data_date, ".png"), width = 800, height = 800)
    # matplot(t(friction_basis$basis_coefs), type = "l")
    plot_domain <- 1:2001 #1000
    plot(friction_basis$true_vals[1, plot_domain], type = "l")
    lines(friction_basis$fitted_values[1, plot_domain], col  = "red")

    dev.off()

    ## De-trend the bedrock
    df <- data.frame(obs_locations = domain[bed_obs$locations], bed_elev = bed_obs$obs)
    bed.fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.25,
                    control = loess.control(surface = "direct"))
    bed_mean <- predict(bed.fit, newdata = data.frame(obs_locations = domain))
    bed_arr_s <- sim_param_list$bedrock

    # bed_mean <- colMeans(bed_arr_s)
    mean_mat <- matrix(rep(bed_mean, nrow(bed_arr_s)), nrow = nrow(bed_arr_s), ncol = length(bed_mean), byrow = T)
    bed_arr_demean <- bed_arr_s - mean_mat

    ## Fit basis to bedrock
    bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, bed_arr = bed_arr_demean,
                                lengthscale = 5e3)
    # bed_fit <- list(mean = bed_mean, basis = bed_basis)
    bed_basis$mean <- bed_mean

    plot_domain <- 1:2001
    plot(bed_mean + bed_basis$true_vals[1, plot_domain], type = "l")
    lines(bed_mean + bed_basis$fitted_values[1, plot_domain], col  = "red")

    if (save_sims) {
      setf <- formatC(set, width = 2, flag = "0")
      qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
      qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
    }

    ## Pair the bed-friction observations into a list
      fitted_param_list <- lapply(1:N, function(r) {
        list(
          # friction = exp(friction_basis$fitted_values[r, ]),
          friction = friction_basis$fitted_values[r, ],
          bedrock = bed_basis$fitted_values[r, ] + bed_mean
        )
      })
    # }

    ## Generate observations based on the simulated bed and friction
    test <- try({
      # sim_results <- run_sims(nsims = N, years = years, sim_beds = T,
      #                         bed_obs = bed_obs, steady_state = ssa_steady)
      sim_results <- sim_obs(
        param_list = fitted_param_list,
        years = years, # sim_beds = T,
        # bed_obs = bed_obs,
        warmup = warmup,
        steady_state = ssa_steady,
        log_transform = log_transform
      )

      sim_results_list[[i]] <- sim_results
    
      })

    if (save_sims) {
      qsave(sim_results, file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))
    }
  }

} else {


  # sim_results_list <- mclapply(sets, function(set) {
  #   setf <- formatC(set, width = 2, flag = "0")
  #   sim_results <- qread(file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))
  # }, mc.cores = 10L)
  #   # sim_results_list[[i]] <- sim_results
    
  # flags <- lapply(sim_results_list, function(sim_results) {

  #   if (length(sim_results$errors) > 0) {
  #     flag <- 1
  #   } else {
  #     flag <- 0
  #   }

  #   return(flag)
  #   }
  # )

}

##########################################
##      Process simulation results      ##
##########################################

if (reprocess_sim_results) {
  for (i in 1:length(sets)) {
    set <- sets[i]

    sim_results <- sim_results_list[[i]]
    
    ## Need to get rid of the simulations that failed here
    bad_sims <- sim_results$bad_sims
    if (length(bad_sims) > 0) {
      cat("Some simulations failed in set", set, "\n")
      good_sims <- sim_results$results[-bad_sims]

      ## Delete corresponding bad simulations from the fitted basis
      friction_basis$basis_coefs <- friction_basis$basis_coefs[-bad_sims, ]
      friction_basis$true_vals <- friction_basis$true_vals[-bad_sims, ]
      friction_basis$fitted_values <- friction_basis$fitted_values[-bad_sims, ]

      bed_basis$basis_coefs <- bed_basis$basis_coefs[-bad_sims, ]
      bed_basis$fitted_values <- bed_basis$fitted_values[-bad_sims, ]
      bed_basis$true_vals <- bed_basis$true_vals[-bad_sims, ]

    } else {
      good_sims <- sim_results$results
    }

    generated_data <- process_sim_results(sims = good_sims)

    surface_obs_arr_s <- generated_data$surface_obs
    friction_arr_s <- generated_data$friction_arr

    true_surface_elevs <- generated_data$true_surface_elevs
    true_thicknesses <- generated_data$true_thicknesses
    true_velocities <- generated_data$true_velocities

    gl_arr_s <- generated_data$gl_arr
    bed_arr_s <- generated_data$bed_arr

    ## Should scale the friction values here
    # friction_arr_s <- friction_arr_s/fric_scale

    # ## Basis function representation of the friction coefficients
    # fitted_friction_s <- fit_basis(nbasis = nbasis, domain = domain, friction_arr = friction_arr_s)

    if (save_sims) {
      setf <- formatC(set, width = 2, flag = "0")
      # qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
      # qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
      qsave(surface_obs_arr_s, file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
      qsave(friction_arr_s, file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
      qsave(gl_arr_s, file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))
      qsave(bed_arr_s, file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
      qsave(true_surface_elevs, file = paste0(train_data_dir, "/true_surface_elevs_", setf, "_", data_date, ".qs"))
      qsave(true_thicknesses, file = paste0(train_data_dir, "/true_thicknesses_", setf, "_", data_date, ".qs"))
      qsave(true_velocities, file = paste0(train_data_dir, "/true_velocities_", setf, "_", data_date, ".qs"))
    }
  }

}



##############################################
##      Plot some simulations to check      ##
##############################################

set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")

### True thickness, friction, bed, grounding line
# thickness_velocity_arr <- readRDS(file = paste0(train_data_dir, "/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
surface_obs_arr <- qread(file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
friction_arr <- qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
bed_arr <- qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))

## Fitted friction and bed
friction_basis <- qread(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
# bed_basis <- readRDS(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".rds"))

bed_basis <- qread(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
# bed_arr <- bed_data$basis$fitted_values
# bed_basis <- bed_data$basis

friction_arr <- friction_basis$true_vals
fitted_friction <- friction_basis$fitted_values
bed_arr <- bed_basis$true_vals
fitted_bed <- bed_basis$fitted_values
bed_mean <- bed_basis$mean

## Plot velocity array for a simulation
png("plots/input/surface_obs_sim.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(2, 1))
matplot(surface_obs_arr[1,,,1], type = "l", main = "Surface elevation (m)", ylab = "Surface elev. (m)", xlab = "Domain (grid points)")
matplot(surface_obs_arr[1,,,2], type = "l", main = "Velocity (m/yr)", ylab = "Velocity (m/yr)", xlab = "Domain (grid points)")
dev.off()

## Hovmoller plots of surface elevation and velocity
plots <- list()

nsamples <- 2
sims <- sample(1:N, size = nsamples)

space <- domain / 1000
time <- 1:((years - 1) - warmup) # 1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples * 4), nsamples, 4, byrow = T)

gl <- ceiling(gl_arr[1, 1] / (domain[length(domain)] / 1000) * length(domain))

# s <- 1
for (s in 1:nsamples) {
  sim <- sims[[s]]
  ind <- inds[s, ]
  # print(ind)
  # thickness_velocity_arr <- sim_results[[s]]$thickness_velocity_arr
  # thickness <- thickness_velocity_arr[sim,,,1]
  # velocity <- thickness_velocity_arr[sim,,,2]
  surface_elev <- surface_obs_arr[sim, , time, 1]
  velocity <- surface_obs_arr[sim, , time, 2]
  # grid_test$thickness <- as.vector(thickness)
  grid_test$surface_elev <- as.vector(surface_elev)
  grid_test$velocity <- as.vector(velocity)

  # thickness_plot <- ggplot(grid_test) +
  # geom_tile(aes(space, time, fill = thickness)) +
  surface_elev_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = surface_elev)) +
    scale_y_reverse(
        breaks = seq(min(grid_test$time), max(grid_test$time), 5),
        labels = seq(min(grid_test$time), max(grid_test$time)+1, 5)
    ) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    xlab("Domain (km)") +
    ylab("Time (yr)") +
    labs(fill = "Surface elev. (m)")

  velocity_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = velocity)) +
    scale_y_reverse(
        breaks = seq(min(grid_test$time), max(grid_test$time), 5),
        labels = seq(min(grid_test$time), max(grid_test$time)+1, 5)
    ) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    xlab("Domain (km)") +
    ylab("Time (yr)") +
    labs(fill = bquote("Velocity (m" ~"yr"^-1 ~ ")"))


  if (log_transform) {
    fitted_fric_sim <- exp(fitted_friction[sim, 1:gl])
    friction_sim <- exp(friction_arr[sim, 1:gl])
  } else {
    fitted_fric_sim <- fitted_friction[sim, 1:gl]
    friction_sim <- friction_arr[sim, 1:gl]
  }

  fric_df <- data.frame(
    domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim,
    fitted_fric = fitted_fric_sim
  )
  friction_plot <- fric_df %>% ggplot() +
    geom_line(aes(x = domain, y = fitted_fric), lwd = 1) +
    # geom_line(aes(x = domain, y = friction), col = "red", lwd = 1) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    xlab("Domain (km)") +
    ylab(bquote("Friction (M Pa m"^"-1/3" ~ "yr"^"1/3" ~ ")"))

  bed_sim <- bed_arr[sim, ] + bed_mean
  fitted_bed_sim <- fitted_bed[sim, ] + bed_mean
  bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim, fitted_bed = fitted_bed_sim)
  
  bed_plot <- bed_df %>% ggplot() + 
    geom_line(aes(x = domain, y = fitted_bed), lwd = 1) +
    # geom_line(aes(x = domain, y = bed), col = "red", lwd = 1) +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    # xlim(0, 800) + 
    # ylim(-1250, -300) +
    xlab("Domain (km)") +
    ylab("Bed (m)")

  # plots[[ind[1]]] <- thickness_plot
  plots[[ind[1]]] <- surface_elev_plot
  plots[[ind[2]]] <- velocity_plot
  plots[[ind[3]]] <- friction_plot
  plots[[ind[4]]] <- bed_plot

}

plots <- lapply(plots, function(p) p + theme(plot.margin = margin(40, 20, 20, 20)))

png(file = paste0("./plots/simulations_", setf, "_", data_date, "_01.png"), 
          width = 4000, height = 2200 * nsamples, res = 300)
# grid.arrange(grobs = plots, ncol = nsamples, nrow = 4)
grid.arrange(grobs = plots, nrow = nsamples, ncol = 4, 
            layout_matrix= matrix(1:(nsamples*4), 4, nsamples))
dev.off()

# ## Plot an example ice sheet profile
# top_surface <- ssa_steady$top_surface
# bottom_surface <- ssa_steady$bottom_surface
# bedrock <- ssa_steady$bedrock

# # timepts <- seq(3001, 2501, -100)
# # top_surf_retreat <- ssa_steady$all_top_surface[, timepts]
# ref_df <- data.frame(domain = ssa_steady$domain / 1000, 
#                     top_surface = top_surface, 
#                     # top_surface_1 = top_surf_retreat[, 1],
#                     # top_surface_2 = top_surf_retreat[, 2],
#                     # top_surface_3 = top_surf_retreat[, 3],
#                     # bottom_surface= bottom_surface, 
#                     bedrock = bedrock)


# ## Plot ice sheet profile
# ss_plot <- ggplot(ref_df, aes(x = domain)) +
#   geom_line(aes(y = top_surface)) +
#   geom_line(aes(y = bottom_surface)) +
#   geom_line(aes(y = bedrock)) +
#   theme_bw() +
#   theme(text = element_text(size = 24)) +
#   xlim(250, 450) +
#   ylim(-1000, 1500) +
#   xlab("Domain (km)") +
#   ylab("Elevation (m)")

# # png(file = paste0("./plots/steady_state_profile.png"), width = 800, height = 500)
# # print(ss_plot)
# # dev.off()

# ## Now plot the evolution of the ice sheet over time
# sim <- 1
# surface_elev <- surface_obs_arr[sim, , , 1]
# # velocity <- surface_obs_arr[sim, , , 2]
# bed <- bed_arr[sim, ]

# # ice_evol_df <- data.frame(domain = ssa_steady$domain / 1000, 
# #                           surface_elev_1 = surface_elev[, 1], 
# #                           surface_elev_2 = surface_elev[, 2],
# #                           surface_elev_3 = surface_elev[, 3],

# #                           bed = bed)

# ref_df$bed_sim <- bed
# ss_plot2 <- ggplot(ref_df, aes(x = domain)) +
#   geom_line(aes(y = top_surface)) +
#   geom_line(aes(y = bottom_surface)) +
#   geom_line(aes(y = bed_sim, col = "red")) +
#   theme_bw() +
#   theme(text = element_text(size = 24)) +
#   xlim(250, 450) +
#   ylim(-1000, 1500) +
#   xlab("Domain (km)") +
#   ylab("Elevation (m)")


# png(file = paste0("./plots/steady_state_profile2.png"), width = 800, height = 500)
# print(ss_plot2)
# dev.off()

# ## Plot fitted frictions
# nsamples <- 2
# par(mfrow = c(nsamples/2, 2))

# for (sim in 1:nsamples) {
#   plot_domain <- 1:1000
#   plot(domain[plot_domain]/1000, friction_arr[sim, plot_domain,,]/fric_scale,
#        type = "l", lwd = 1.5, xlab = "Domain (km)", ylab = "Friction (unit)")
#   # lines(domain[plot_domain]/1000, lmfit$fitted.values[plot_domain]/fric_scale, col = "seagreen", lwd = 1.5)

#   lines(domain[plot_domain]/1000, fitted_friction[sim, plot_domain], col = "red", lwd = 1.5)
#   # legend("topright", legend = c("global basis", "local basis"), col = c("seagreen", "red"), lty = 1, lwd = 1.5)
#   # legend("topright", legend = c("original friction", "local basis rep"), col = c("black", "red"), lty = 1, lwd = 1.5)

# }



