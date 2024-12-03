## Generate 1000 simulations

setwd("~/SSA_model/CNN/real_data/")

rm(list = ls())

library("parallel")
library("Matrix")
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
library("R.utils")
# library("sp")
library("fields")
# library("tidyr")
library("dplyr")
library("matrixStats") # for the rowMaxs() function
library("mvtnorm")
library(abind)
library("ggplot2")
# library("plotly")
library("gridExtra")
library("FRK")
library(qs)

source("./source/sim_params.R")
# source("./source/run_sims.R")
source("./source/sim_obs.R")
source("./source/cond_sim_gp.R")
source("./source/get_ini_thickness.R")
source("./source/process_sim_results.R")
source("./source/fit_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
# source("./source/create_ref.R")
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/get_surface_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")


## Some flagsn,
regenerate_sims <- T
refit_basis <- T
save_sims <- T
log_transform <- T
# sim_beds <- Tv

train_data_dir <- "./data/training_data"

## Presets
data_date <- "20241111" #"20241103" 
N <- 1000 # number of simulations per set
sets <- 1#:50 #:10
setf <- paste0("sets", sets[1], "-", sets[length(sets)])

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 70 #150

# 0. Load ice sheet at steady state
ssa_steady <- qread(file = paste0(train_data_dir, "/steady_state/steady_state_", data_date, ".qs"))
surf_elev_mat <- qread(file = "./data/surface_elev/surf_elev_mat.qs")

## 0. Flowline
flowline <- readRDS("./data/flowline_regrid.rds")
# J <- 2001
# flowline <- flowline[1:J, ]
J <- nrow(flowline)
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
domain <- c(0, cumsum(na.omit(flowline_dist)))

# 0. Read bed observations
bed_obs_df <- qread(file = paste0("./data/bed_obs_df.qs"))

years <- 20
domain <- ssa_steady$domain

## Scaling units for friction coefficients
secpera <- 31556926
fric_scale <- 1e6 * secpera^(1 / 3)

t1 <- proc.time()
if (regenerate_sims) {

  flags <- c()

  bed_sim_list <- list()
  fric_sim_list <- list()

  for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Generating set", set, "\n")
    
    sim_param_output <- sim_params(
      nsims = N, domain = ssa_steady$domain,
      bed_obs = bed_obs_df[bed_obs_df$chosen == 1, ]
    )

    # print("Simulating bed elevations...")
    # bed_sims <- simulate_bed(N, domain = domain, 
    #                         obs_locations = bed_obs$ind, 
    #                         obs = bed_obs$bed_elev) 
    sim_param_list <- sim_param_output$sim_param_list
    bed_sims <- sim_param_list$bedrock
    fric_sims <- sim_param_list$friction

    bed_sim_list[[i]] <- bed_sims
    fric_sim_list[[i]] <- fric_sims

    ## Plot conditional bed simulations
    bed_sims_df <- data.frame(domain = domain, 
                              bed1 = bed_sims[1, ], 
                              bed2 = bed_sims[2, ],
                              bed3 = bed_sims[3, ],
                              bed4 = bed_sims[4, ],
                              bed5 = bed_sims[5, ],
                              bed_obs = bed_obs_df$bed_elev, 
                              obs_loc = bed_obs_df$loc) #,
                              # chosen = bed_obs_df$chosen)

    png(file = paste0("./plots/bed/bed_sims_", setf, "_", data_date, ".png"), width = 1000, height = 500)
    bed_sim_plot <- bed_sims_df %>% ggplot() + 
      geom_line(aes(x = domain, y = bed1), col = "grey") +
      geom_line(aes(x = domain, y = bed2), col = "grey") +
      geom_line(aes(x = domain, y = bed3), col = "grey") +
      geom_line(aes(x = domain, y = bed4), col = "grey") +
      geom_line(aes(x = domain, y = bed5), col = "grey") +
      geom_point(data = bed_obs_df, aes(x = loc, y = bed_elev)) +
      theme_bw()
    print(bed_sim_plot)
    dev.off()

    if (save_sims) {
      qsave(bed_sims, file = paste0(train_data_dir, "/bed_sims_", data_date, ".qs"))
      qsave(fric_sims, file = paste0(train_data_dir, "/fric_sims_", data_date, ".qs"))
    }
  }

  ## Concatenate the bed simulations and take the mean
  bed_sims_arr <- abind(bed_sim_list, along = 1)
  # fric_sims_arr <- abind(fric_sim_list, along = 1)
  bed_mean <- colMeans(bed_sims_arr)
    
  ## Then fit basis here  
  for (i in 1:length(sets)) {
    set <- sets[i]
    cat("Fitting basis functions for set", set, "\n")
    setf <- formatC(set, width = 2, flag = "0")

    bed_sims <- bed_sim_list[[i]] 
    fric_sims <- fric_sim_list[[i]] 

    ## Fit basis to log(friction)
    friction_basis <- fit_friction_basis(
      nbasis = nbasis,
      domain = domain,
      fric_arr = fric_sims,
      log_transform = log_transform#,
      # lengthscale = 10e3
    )

    ## De-trend the bedrock
    # bed_mean <- colMeans(bed_sims)
    mean_mat <- matrix(rep(bed_mean, nrow(bed_sims)), nrow = nrow(bed_sims), ncol = length(bed_mean), byrow = T)
    bed_arr_demean <- bed_sims - mean_mat

    ## Fit basis to bedrock
    bed_basis <- fit_bed_basis(nbasis = nbasis, domain = domain, 
                              bed_arr = bed_arr_demean) #,
                              # lengthscale = 10e3)
    # bed_fit <- list(mean = bed_mean, basis = bed_basis)
    bed_basis$mean <- bed_mean

    png(file = paste0("./plots/friction/friction_basis_", setf, "_", data_date, ".png"))
    plot_domain <- 1:J #1000
    plot(domain[plot_domain], friction_basis$true_vals[2, plot_domain], type = "l")
    lines(domain[plot_domain], friction_basis$fitted_values[2, plot_domain], col  = "red")
    dev.off()

    png(file = paste0("./plots/bed/bed_basis_", setf, "_", data_date, ".png"))
    plot_domain <- 1:J #1000
    plot(domain[plot_domain], bed_mean + bed_basis$true_vals[2, plot_domain], type = "l")
    # matplot(domain[plot_domain], t(bed_sims), type = "l")
    plot(domain[plot_domain], bed_mean + bed_sims[2, ], type = "l")
    # lines(domain[plot_domain], bed_sims[5, ], col  = "red")
    # lines(domain[plot_domain], bed_mean + bed_basis$fitted_values[8, plot_domain], col  = "red")
    dev.off()

    ## Pair the bed-friction observations into a list
    fitted_param_list <- lapply(1:N, function(r) {
      list(
        # friction = exp(friction_basis$fitted_values[r, ]),
        friction = friction_basis$fitted_values[r, ],
        bedrock = bed_basis$fitted_values[r, ] + bed_mean
      )
    })

    # 0. Velocity and ice thickness data
    ## Initial thickness
    # se_grounded <- na.omit(surf_elev_mat[, 21]) # Use surface elevation at the final time to initialise ice thickness
    # H_ini <- se_grounded - bed_sims[, 1:length(se_grounded)]
    # missing <- which(is.na(surf_elev_mat[, 21]))
    # rho <- 910.0
    # rho_w <- 1028.0
    # thickness_at_tail <- - bed_sims[missing, ] * rho_w / rho # worked out based on grounding line conditions
    # H_ini_new <- c(H_ini, thickness_at_tail) #+ offset
    # H_ini_new <- H_ini_new + offset

    ini_thickness <- #lapply(1:N, function(r) {
      get_ini_thickness(surf_elev = surf_elev_mat[, 21], bed = colMeans(bed_sims))
    # })

# png("./plots/temp/ini_thickness.png")
# plot(domain, ini_thickness, type = "l")
# dev.off()

    ## Velocity
    vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
    vel_curr <- vel_mat[, ncol(vel_mat)]
    ## Smooth the velocity out with loess
    vel_curr_smooth <- loess(vel_curr ~ domain, span = 0.1)

    ## Generate observations based on the simulated bed and friction
    test <- try(
      # sim_results <- run_sims(nsims = N, years = years, sim_beds = T,
      #                         bed_obs = bed_obs, steady_state = ssa_steady)
      sim_results <- sim_obs(
        param_list = fitted_param_list,
        years = years, # sim_beds = T,
        # bed_obs = bed_obs,
        # steady_state = ssa_steady,
        ini_thickness = ini_thickness,
        ini_velocity = vel_curr_smooth$fitted,
        log_transform = log_transform
      )
    )

    # calc <- if(length(sim_results$errors) > 0) {
    #   flags[i] <- 1
    #   next
    # } else {
    #   flags[i] <- 0
    # }

    if (save_sims) {
      qsave(sim_results, file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))
    }

    # sim_results <- readRDS(file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".rds"))

    ## Need to get rid of the simulations that failed here
    bad_sims <- sim_results$bad_sims
    # sim_results$errors[[1]]
    # ind <- bad_sims[1]
    # bad_bed <- sim_results$params[[ind]]$bedrock
    # bad_fric <- sim_results$params[[ind]]$friction
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

    # thickness_velocity_arr_s <- generated_data$thickness_velocity_arr
    surface_obs_arr_s <- generated_data$surface_obs_arr
    friction_arr_s <- generated_data$friction_arr
    gl_arr_s <- generated_data$gl_arr
    bed_arr_s <- generated_data$bed_arr

    ## Should scale the friction values here
    # friction_arr_s <- friction_arr_s/fric_scale

    true_surface_elevs <- generated_data$true_surface_elevs
    true_thicknesses <- generated_data$true_thicknesses
    true_velocities <- generated_data$true_velocities


    if (save_sims) {
      setf <- formatC(set, width = 2, flag = "0")
      qsave(friction_basis, file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))
      qsave(bed_basis, file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
      # qsave(thickness_velocity_arr_s, file = paste0(train_data_dir, "/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
      qsave(surface_obs_arr_s, file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))
      qsave(friction_arr_s, file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
      qsave(gl_arr_s, file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))
      qsave(bed_arr_s, file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
      qsave(true_surface_elevs, file = paste0(train_data_dir, "/true_surface_elevs_", setf, "_", data_date, ".qs"))
      qsave(true_thicknesses, file = paste0(train_data_dir, "/true_thicknesses_", setf, "_", data_date, ".qs"))
      qsave(true_velocities, file = paste0(train_data_dir, "/true_velocities_", setf, "_", data_date, ".qs"))
    }
  }
} else {

  # flags <- c()
  # for (s in 1:length(sets)) {
  #   set <- sets[s]
  #   cat("Checking set", set, "\n")
  #   setf <- formatC(set, width = 2, flag = "0")
  #   sim_results <- qread(file = paste0(train_data_dir, "/sim_results_", setf, "_", data_date, ".qs"))

  #   if (length(sim_results$errors) > 0) {
  #     flags[s] <- 1
  #   } else {
  #     flags[s] <- 0
  #   }
  # }

# bad_sets <- sets[which(flags == 1)]
# sink("bad_sets.txt")
# print(bad_sets)
# sink()
}


t2 <- proc.time()

## Check if the generated data is close to real data here
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs")

png("./plots/temp/surf_elev.png")
matplot(surface_obs_arr[1,,,1], type = "l", col = "grey")
matlines(surf_elev_mat, col = "red", lty = 2)
dev.off()

browser()
## Plot some simulations to check
set <- sets[1]
setf <- formatC(set, width = 2, flag = "0")

### True thickness, friction, bed, grounding line
# thickness_velocity_arr <- readRDS(file = paste0(train_data_dir, "/thickness_velocity_arr_", setf, "_", data_date, ".rds"))
surface_obs_arr <- qread(file = paste0(train_data_dir, "/surface_obs_arr_", setf, "_", data_date, ".qs"))

# friction_arr <- qread(file = paste0(train_data_dir, "/friction_arr_", setf, "_", data_date, ".qs"))
# bed_arr <- qread(file = paste0(train_data_dir, "/bed_arr_", setf, "_", data_date, ".qs"))
gl_arr <- qread(file = paste0(train_data_dir, "/gl_arr_", setf, "_", data_date, ".qs"))

## Fitted friction and bed
friction_basis <- qread(file = paste0(train_data_dir, "/friction_basis_", setf, "_", data_date, ".qs"))

# bed_basis <- readRDS(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".rds"))

bed_basis <- qread(file = paste0(train_data_dir, "/bed_basis_", setf, "_", data_date, ".qs"))
# bed_arr <- bed_data$basis$fitted_values
# bed_basis <- bed_data$basis

friction_arr <- friction_basis$true_vals
# bed_mean_mat <- matrix(rep(bed_basis$mean, nrow(bed_basis$true_vals)), nrow = nrow(bed_basis$true_vals), ncol = length(bed_basis$mean), byrow = T)
bed_arr <- bed_basis$true_vals #+ bed_mean_mat
fitted_friction <- friction_basis$fitted_values
fitted_bed <- bed_basis$fitted_values
bed_mean <- bed_basis$mean

# fitted_friction <- fitted_friction[-bad_sims, ]
# fitted_bed <- fitted_bed[-bad_sims, ]

plots <- list()

nsamples <- 4
sims <- sample(1:dim(surface_obs_arr)[1], size = nsamples)

space <- domain / 1000
time <- 1:(years + 1) # 1:dim(thickness_velocity_arr)[3]
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
  surface_elev <- surface_obs_arr[sim, , , 1]
  velocity <- surface_obs_arr[sim, , , 2]
  # grid_test$thickness <- as.vector(thickness)
  grid_test$surface_elev <- as.vector(surface_elev)
  grid_test$velocity <- as.vector(velocity)

  # thickness_plot <- ggplot(grid_test) +
  # geom_tile(aes(space, time, fill = thickness)) +
  surface_elev_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = surface_elev)) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    theme_bw() +
    # labs(fill="Thickness (m)")
    labs(fill = "Surface elevation (m)")

  velocity_plot <- ggplot(grid_test) +
    geom_tile(aes(space, time, fill = velocity)) +
    theme_bw() +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    labs(fill = bquote("Velocity (m" ~ a^-1 ~ ")"))


  if (log_transform) {
    fitted_fric_sim <- exp(fitted_friction[sim, 1:gl])
    friction_sim <- exp(friction_arr[sim, 1:gl])
  } else {
    fitted_fric_sim <- fitted_friction[sim, 1:gl]
    friction_sim <- friction_arr[sim, 1:gl]
  }

  df <- data.frame(
    domain = ssa_steady$domain[1:gl] / 1000, friction = friction_sim,
    fitted_fric = fitted_fric_sim
  )
  friction_plot <- ggplot(df, aes(x = domain, y = friction)) +
    geom_line() +
    geom_line(aes(x = domain, y = fitted_fric), col = "red") +
    theme_bw() +
    xlab("Domain (km)") +
    ylab(bquote("Friction (M Pa m"^"-1/3" ~ "a"^"1/3" ~ ")"))

  bed_sim <- bed_arr[sim, ] + bed_mean
  fitted_bed_sim <- fitted_bed[sim, ] + bed_mean
  bed_df <- data.frame(domain = ssa_steady$domain / 1000, bed = bed_sim, fitted_bed = fitted_bed_sim)
  
  bed_plot <- ggplot(bed_df, aes(x = domain, y = bed)) +
    geom_line() +
    geom_line(aes(x = domain, y = fitted_bed), col = "red") +
    theme_bw() +
    xlab("Domain (km)") +
    ylab("Bed (m)")

  # plots[[ind[1]]] <- thickness_plot
  plots[[ind[1]]] <- surface_elev_plot
  plots[[ind[2]]] <- velocity_plot
  plots[[ind[3]]] <- friction_plot
  plots[[ind[4]]] <- bed_plot

}

png(file = paste0("./plots/simulations_", setf, "_", data_date, ".png"), width = 2000, height = 800)
grid.arrange(grobs = plots, nrow = nsamples, ncol = 4)
dev.off()


png(file = paste0("./plots/cnn/test.png"))
matplot(surface_obs_arr[2,,,1], type = "l", col = "grey")
matlines(surf_elev_mat, col = "red")
dev.off()


