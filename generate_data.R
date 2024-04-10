## Generate 1000 simulations 

setwd("/home/babv971/SSA_model/solveSSA")

rm(list = ls())

library("parallel")
library("Matrix")
# library("qlcMatrix")
# library("fastmatrix")
# library("expm")
# library("R.utils")
# library("sp")
library("fields")
# library("tidyr")
library("dplyr")
library("matrixStats") # for the rowMaxs() function
library("mvtnorm")
library("parallel")

library("ggplot2")
# library("plotly")
library("gridExtra")
library(FRK)

# library(gstat)
# library(sp)
# library("mvnfast")
# library("splines")

source("./source/run_sims.R")
source("./source/fit_fric_basis.R")
source("./source/surface_elev.R")
source("./source/create_params.R")
source("./source/create_ref.R")
source("./source/solve_ssa_nl.R")
# source("./source/solve_ssa_nl2.R")

source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
source("./source/mvnorm_sample.R")
source("./source/get_obs.R")
# source("simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")
# source("initialise_ens.R")
# source("ssa_plot_ini_ens.R")
source("./source/propagate.R")
source("./source/obs_operator.R")
# source("run_enkf.R")
# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
# source("initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# ## Seed for generating bed
# ssa_seed <- 123
# set.seed(ssa_seed)

## Some flags
regenerate_sims <- T
save_sims <- T

## Presets
data_date <- "20240320" #"20220329" 
N <- 10000 # number of simulations per set
sets <- 6:10

# set <- 1 #commandArgs(trailingOnly = TRUE)

nbasis <- 150

# 0. Load ice sheet at steady state
ssa_steady <- readRDS(file = paste("./training_data/initial_conds/ssa_steady_20220329.rds", sep = ""))

years <- 50
domain <- ssa_steady$domain

secpera <- 31556926
fric_scale <- 1e6 * secpera^(1/3)

if (regenerate_sims) {

  for (set in sets) {
    cat("Generating set", set, "\n")
    try(generated_data <- run_sims(nsims = N)) 

    thickness_velocity_arr_s <- generated_data$thickness_velocity_arr
    friction_arr_s <- generated_data$friction_arr
    gl_arr_s <- generated_data$gl_arr

    fitted_friction_s <- fit_fric_basis(nbasis = nbasis, domain = domain, friction_arr = friction_arr_s)

    if (save_sims) {
        saveRDS(thickness_velocity_arr_s, file = paste0("./training_data/thickness_velocity_arr_", set, "_", data_date))
        saveRDS(friction_arr_s, file = paste0("./training_data/friction_arr_", set, "_", data_date))
        saveRDS(gl_arr_s, file = paste0("./training_data/gl_arr_", set, "_", data_date))
        saveRDS(fitted_friction_s, file = paste0("./training_data/fitted_friction_", set, "_", data_date))
    }
  }
  
} 


## Post processing
set <- 1
thickness_velocity_arr <- readRDS(file = paste0("./training_data/thickness_velocity_arr_", set, "_", data_date))
friction_arr <- readRDS(file = paste0("./training_data/friction_arr_", set, "_", data_date))
gl_arr <- readRDS(file = paste0("./training_data/gl_arr_", set, "_", data_date))
fitted_friction <- readRDS(file = paste0("./training_data/fitted_friction_", set, "_", data_date))
# if (regenerate_sims) {

# # 1. Simulate friction coefficients
#   print("Simulating friction coefficient...")
#   # simulated_friction <- simulate_friction(nsim = N, domain = ssa_steady$domain,
#                                           # sill = 8e-5, nugget = 0, range = 2.5e3) # doesn't need input, or maybe just put the mean in there
  
#   fric.sill <- 8e-5
#   fric.nugget <- 0
#   fric.range <- 10e3
  
#   simulated_friction <- simulate_friction2(nsim = N, domain = domain,
#                                           sill = fric.sill, nugget = fric.nugget, 
#                                           range = fric.range) # doesn't need input, or maybe just put the mean in there
#   sim_fric_list <- lapply(1:N, function(c) simulated_friction[, c])
#   # 1. Reference run
#   print("Simulating observations...")

#   # sim_results <- list()
#   grounding_lines <- list()
  
#   ## Process noise parameters
#   ones <- rep(1, length(domain))
#   D <- rdist(domain)
#   l <- 50e3
#   R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
#   L <- t(chol(R))
#   L <- as(L, "dgCMatrix")
#   process_noise_info <- list(corrmat_chol = L, length_scale = l)

#   t1 <- proc.time()
#   sim_results <- mclapply(sim_fric_list, function(simulated_friction) {
#     # cat("sim =", sim, "\n")
#     # steady_state = ssa_steady
#     reference <- solve_ssa_nl(domain = ssa_steady$domain, 
#                               bedrock = ssa_steady$bedrock, 
#                               friction_coef = simulated_friction,
#                               ini_velocity = ssa_steady$current_velocity,
#                               ini_thickness = ssa_steady$current_thickness,
#                               years = years, steps_per_yr = 52, 
#                               save_model_output = TRUE, 
#                               perturb_hardness = TRUE,
#                               add_process_noise = T,
#                               process_noise_info = process_noise_info)    
    
#     thickness_velocity_obs <- array(data = cbind(reference$all_thicknesses[, 2:(years+1)], 
#                                 reference$all_velocities[, 2:(years+1)]),
#                   dim = c(length(domain), years, 2))
#     gl <- reference$grounding_line
    
#     simulated_data <- list(thickness_velocity_arr = thickness_velocity_obs,
#                            friction_arr = simulated_friction,
#                            grounding_line = gl)
#     # simulated_data <- obs
#     return(simulated_data)
#   },
#   #steady_state = ssa_steady, 
#   mc.cores = 50L)
#   t2 <- proc.time()


#  ## Note to self: should save a version with the "true" friction coef as well
#   thickness_velocity_list <- lapply(1:N, function(i) sim_results[[i]]$thickness_velocity_arr)
#   friction_list <- lapply(1:N, function(i) sim_results[[i]]$friction_arr)
#   gl_list <- lapply(1:N, function(i) sim_results[[i]]$grounding_line)

#   concat_input <- do.call("rbind", thickness_velocity_list)
#   thickness_velocity_arr <- array(concat_input, dim = c(N, dim(thickness_velocity_list[[1]])))

#   concat_output <- do.call("rbind", friction_list)
#   friction_arr <- array(concat_output, dim = c(N, length(friction_list[[1]]), 1L, 1L))

#   concat_gl <- do.call("rbind", gl_list)
#   gl_arr <- array(concat_gl, dim = c(N, 1L, length(gl_list[[1]]), 1L))

#   if (save_sims) {
#     saveRDS(thickness_velocity_arr, file = paste0("./output/thickness_velocity_arr_", set, "_", data_date))
#     saveRDS(friction_arr, file = paste0("./output/friction_arr_", set, "_", data_date))
#     saveRDS(gl_arr, file = paste0("./output/gl_arr_", set, "_", data_date))
#   }
# } else {
  # thickness_velocity_arr <- readRDS(file = paste0("./output/thickness_velocity_arr_", set, "_", data_date))
  # friction_arr <- readRDS(file = paste0("./output/friction_arr_", set, "_", data_date))
  # gl_arr <- readRDS(file = paste0("./output/gl_arr_", set, "_", data_date))
# }


# ## Plots
plots <- list()

nsamples <- 2
sims <- sample(1:N, size = nsamples)

space <- domain/1000
time <- 1:years #1:dim(thickness_velocity_arr)[3]
grid_test <- expand.grid(space, time)
head(grid_test)
names(grid_test) <- c("space", "time")

inds <- matrix(1:(nsamples*3), nsamples, 3, byrow = T)

gl <- ceiling(gl_arr[1,,1,]/(domain[length(domain)]/1000)*length(domain))

for (s in 1:nsamples) {
  sim <- sims[[s]]
  ind <- inds[s, ]
  # print(ind)
  # thickness_velocity_arr <- sim_results[[s]]$thickness_velocity_arr
  thickness <- thickness_velocity_arr[sim,,,1]
  velocity <- thickness_velocity_arr[sim,,,2]
  grid_test$thickness <- as.vector(thickness)
  grid_test$velocity <- as.vector(velocity)

  thickness_plot <- ggplot(grid_test) + 
    geom_tile(aes(space, time, fill = thickness)) +
    scale_fill_distiller(palette = "Blues", direction = 1) + 
    theme_bw() +
    labs(fill="Thickness (m)")
  velocity_plot <- ggplot(grid_test) + 
    geom_tile(aes(space, time, fill = velocity)) +
    theme_bw() +
    scale_fill_distiller(palette = "Reds", direction = 1) + 
    labs(fill=bquote('Velocity (m'~a^-1~')'))
  
  friction_sim <- friction_arr[sim, 1:gl, , ]
  df <- data.frame(domain = ssa_steady$domain[1:gl]/1000, friction = friction_sim/fric_scale)
  friction_plot <- ggplot(df, aes(x = domain, y = friction)) + geom_line() + 
    theme_bw() + xlab("Domain (km)") + ylab(bquote('Friction (M Pa m'^'-1/3'~'a'^'1/3'~')'))
  
  plots[[ind[1]]] <- thickness_plot
  plots[[ind[2]]] <- velocity_plot
  plots[[ind[3]]] <- friction_plot
  
}

grid.arrange(grobs = plots, nrow = nsamples, ncol = 3)
# ggsave(file = "surfaceplot_2p5e3.png",
#        arrangeGrob(grobs = plots, nrow = nsamples, ncol = 3),
#        width = 4000, height = 1000, units = "px",
#        # width = 4000, height = 3000, units = "px",
#        path = "./output/")

# ## Represent the friction coefficients in terms of basis functions
# nbasis <- 150
# fric.sigma <- fric_cov2(sqrt(fric.sill), fric.range, si = domain, sj = domain)
# eig <- eigen(fric.sigma)
# eofs <- eig$vectors[, 1:nbasis]
# 
# matplot(domain/1000, eofs[, 1:10], type = "l", lwd = 1.5, xlab = "Domain(km)", ylab = "EOFs")
# 
# testsim <- 1
# 
# A <- eig$vectors#eofs
# b <- simulated_friction[, testsim]
# x <- solve(A, b)
# 
# basis50 <- A[, 1:50] %*% x[1:50]
# basis75 <- A[, 1:75] %*% x[1:75]
# basis100 <- A[, 1:100] %*% x[1:100]
# 
# par(mfrow = c(1,1))
# # plot_domain <- 1:1000
# # plot(domain[plot_domain]/1000, simulated_friction[plot_domain, testsim]/fric_scale, 
# #      type = "l", lwd = 1.5,
# #      xlab = "Domain (km)", ylab = "Friction (unit)")
# # lines(domain[plot_domain]/1000, basis50[plot_domain]/fric_scale, col = "red", lwd = 1.5)
# # # lines(domain[plot_domain]/1000, basis75[plot_domain]/fric_scale, col = "purple") # 100 seems better than 50
# # lines(domain[plot_domain]/1000, basis100[plot_domain]/fric_scale, col = "royalblue", lwd = 1.5) # 100 seems better than 50
# # legend("topright", legend = c("K = 50", "K = 100"), 
# #        col = c("red", "royalblue"), lty = 1, lwd = 1.5)
# 
# 
# ## Use lm() to fit a model instead
# df <- as.data.frame(cbind(simulated_friction[, testsim], eofs[, 1:nbasis]))
# colnames(df) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
# lmfit <- lm(fric ~ ., data = df)
# 
# # plot_domain <- 1:1000
# # plot(domain[plot_domain]/1000, simulated_friction[plot_domain, testsim]/fric_scale, 
# #      type = "l", lwd = 1.5, xlab = "Domain (km)", ylab = "Friction (unit)")
# # lines(domain[plot_domain]/1000, lmfit$fitted.values[plot_domain]/fric_scale, col = "seagreen", lwd = 1.5)
# # lines(domain[plot_domain]/1000, basis50[plot_domain]/fric_scale, col = "red", lty = 2, lwd = 1.5)


# ## FRK
# library(FRK)

# nbasis <- 150
# basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis+2)
# basis_centres <- basis_centres[2:(length(basis_centres)-1)]
                        
# testbasis <- local_basis(manifold = real_line(), 
#                           loc = matrix(basis_centres),
#                           type = "bisquare",
#                           scale = rep(5e3, nbasis))
# show_basis(testbasis)
# basis_fns <- lapply(testbasis@fn, function(f) f(domain))
# # plot(domain, basis_fns[[10]])
# basismat <- do.call(cbind, basis_fns)
# matplot(domain/1000, basismat, type = "l", col= "salmon", lty = 1, lwd = 1.5,
#         xlab = "Domain (km)", xlim = c(0, 200))


# fitted_fric <- matrix(NA, length(domain), N)
# for (sim in 1:N) {
#   df_local <- as.data.frame(cbind(simulated_friction[, sim], basismat))
#   colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
#   lmfit_local <- lm(fric ~ ., data = df_local)
#   fitted_fric[, sim] <- lmfit_local$fitted.values
# }

# if (save_sims) {
#   saveRDS(fitted_fric, file = paste0("./output/fitted_fric_", set, "_", data_date))
# }

## Plot fitted frictions
nsamples <- 10
par(mfrow = c(nsamples/2, 2))

for (sim in 1:nsamples) {
  plot_domain <- 1:1000
  plot(domain[plot_domain]/1000, friction_arr[sim, plot_domain,,]/fric_scale,
       type = "l", lwd = 1.5, xlab = "Domain (km)", ylab = "Friction (unit)")
  # lines(domain[plot_domain]/1000, lmfit$fitted.values[plot_domain]/fric_scale, col = "seagreen", lwd = 1.5)
  
  lines(domain[plot_domain]/1000, fitted_friction[sim, plot_domain]/fric_scale, col = "red", lwd = 1.5)
  # legend("topright", legend = c("global basis", "local basis"), col = c("seagreen", "red"), lty = 1, lwd = 1.5)
  # legend("topright", legend = c("original friction", "local basis rep"), col = c("black", "red"), lty = 1, lwd = 1.5)
  
}

