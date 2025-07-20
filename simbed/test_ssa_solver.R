### Main file ###

setwd("/home/babv971/SSA_model/CNN/simbed/")

rm(list = ls())

# library("parallel")
library("Matrix")
library("qlcMatrix")
library("fastmatrix")
library("expm")
library("R.utils")
# library("sp")
# library("fields")
# library("tidyr")
# library("dplyr")
# library("ggplot2")
library("matrixStats") # for the rowMaxs() function
library("mvtnorm")
# library("mvnfast")
# library("splines")
# library("fda")

source("./source/surface_elev.R")
source("./source/create_params.R")
source("./source/create_steady_state.R")
source("./source/create_ref.R")
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")
# source("mvnorm_sample.R")
# source("get_obs.R")
source("./source/simulate_bed.R")
source("./source/simulate_friction.R")
source("./source/azm_cond_sim.R")
# source("initialise_ens.R")
# source("ssa_plot_ini_ens.R")
# source("propagate.R")
# source("obs_operator.R")
# source("run_enkf.R")
# source("run_bg_ens.R")
# source("run_pf.R")
# source("construct_bed_basis.R")
# source("initialise_ice_thickness.R")
# source("compute_block_weights.R")
# source("create_pf_taper.R")
# source("create_pf_smoother.R")

# source("ssa_enkf_plots.R")

## Seed for generating bed
ssa_seed <- 123
set.seed(ssa_seed)

## Some flags
# add_process_noise <- T
# avg_prev_velocity <- TRUE # use the ensemble mean previous velocity to compute the current velocity
# use_ref_velocity <- F # use reference velocity as the initial previous velocity
# use_velocity_dependent_R <- FALSE

## Presets
data_date <- "20220329" #"20230518"
output_date <- "20250402" #"20240518"
rewrite_steady_state <- F
# Ne <- 50 # Ensemble size
# years <- 40
# steps_per_yr <- 100
# n_params <- 20 #number of beds
# n_bed_obs <- 100
# smoothing_factor <- 0.5 # for the PF

# 1. Reference run
print("Simulating ground truth...")

## Physical params
params <- list(
    secpera = 31556926, #seconds per annum
    n = 3.0, # exponent in Glen's flow law
    rho_i = 917.0, #ice density
    rho_w = 1028.0, #sea water density
    #r <- rho / rho_w # define this ratio for convenience
    g = 9.81, #gravity constant
    A = 4.227e-25, #1.4579e-25, # flow rate parameter
    as = 0.28, # surface accumulation rate,
    ab = 0 # basal melt rate, 
)

params$B <- params$A^(-1/params$n)
params$m <- 1/params$n

## Ice shelf specifications
L <- 200e3 # length of domain

## Boundary conditions at the ice divide (x = 0)
# x0 <- 0
# u0 <- 100 / params$secpera # (m/s)
# H0 <- 2000 # (m)

# Discretise domain
J <- 1000 # number of grid 'spaces'
dx <- L / J # grid spacing
x <- seq(0, L, length.out = J+1) # grid point locations

### Friction
# C <- create_fric_coef(x, L) * 1e6 * (params$secpera)^params$m
C <- 7e6

### Bed
bed_params <- list(sill_min = 2000e3, sill_max = 2100e3, sill_slope = 1e-3, b0 = -100, bx = -1e-3)

create_linear_bed <- function(x, params) { # from Alex Robel
  xsill <- x > params$sill_min & x < params$sill_max
  xdsill <- x >= params$sill_max
  sill_length <- params$sill_max - params$sill_min
  
  b <- params$b0 + params$bx * x
  
  b[xsill] <- params$b0 + (params$bx * params$sill_min) + params$sill_slope * (x[xsill] - params$sill_min)
  
  b[xdsill] <- params$b0 + (params$bx * params$sill_min) + params$sill_slope * sill_length +
               params$bx * (x[xdsill] - params$sill_max)
  
  return(b)
}

# b <- create_bed(x) # my own function
b <- create_linear_bed(x, bed_params)

## Initial thickness and velocity

# Grid parameters (from Alex Robel's code)
params$tfinal <- 10e3 * params$year  # Length of transient simulation
params$Nt <- 1e2                     # Number of time steps
params$dt <- params$tfinal / params$Nt # Time step length

params$Nx <- J+1                    # Number of grid points
params$N1 <- floor(J * 0.97)                     # Number of grid points in coarse domain

params$sigGZ <- 0.97                 # Extent of coarse grid (where GL is at sigma=1)

sigma1 <- seq(params$sigGZ / (params$N1 + 0.5), params$sigGZ, length.out = params$N1)
sigma2 <- seq(params$sigGZ, 1, length.out = params$Nx - params$N1 + 1)

params$sigma <- c(sigma1, sigma2[-1])  # Grid points on velocity (includes GL, not ice divide)
params$sigma_elem <- c(0, (params$sigma[1:(params$Nx-1)] + params$sigma[2:params$Nx]) / 2) # Grid points on thickness (includes divide, not GL)
params$dsigma <- diff(params$sigma)    # Grid spacing

params$hscale <- 1000
params$ascale <- 0.1/params$secpera
params$uscale <- (params$rho_i * params$g * params$hscale * params$ascale / C)^(1/(params$m+1)) # velocity scaling
params$xscale <- params$uscale * params$hscale / params$ascale #horizontal distance scaling
params$lambda = 1 - (params$rho_i / params$rho_w) #density difference (lambda param Schoof 2007)
params$sigma <- c(sigma1, sigma2[-1])
params$sigma_elem <- c(0, (params$sigma[1:(params$Nx-1)] + params$sigma[2:params$Nx]) / 2)

params$as <- 1/params$secpera
xg <- 200e3/params$xscale
hf <- (-create_linear_bed(xg * params$xscale, bed_params)/params$hscale)/(1-params$lambda)
h <- 1 - (1 - hf) * params$sigma # initial thickness
u <- 0.3 * (params$sigma_elem^(1/3)) + 1e-3 # initial velocity

h <- h * params$hscale
u <- u * params$uscale * params$secpera

## Now I need to artificially add an ice shelf
shelf_domain <- seq(x[length(x)] + dx, x[length(x)] + 100e3, by = dx)
domain <- c(x, shelf_domain)

## Extend the bedrock
b_shelf <- create_linear_bed(shelf_domain, bed_params)
b <- c(b, b_shelf)

## Extend friction
C_shelf <- 0
friction_coef <- c(rep(C, J), rep(C_shelf, length(shelf_domain)))

## Extend the velocity and thickness
h_shelf <- rep(h[length(h)], length(shelf_domain))
u_shelf <- rep(u[length(u)], length(shelf_domain))

h <- c(h, h_shelf)
u <- c(u, u_shelf)

## Extend the thickness and velocity
ssa_steady <- create_steady_state(use_stored_steady_state = F,
                        add_process_noise_in_ref = F, 
                        rewrite_steady_state = F,
                        # years = 950,
                        phys_params = params,
                        domain = domain, 
                        ini_thickness = h,
                        ini_velocity = u,
                        friction_coef = friction_coef,
                        bedrock = b) 

## Extract individual components
steady_velocity <- ssa_steady$current_velocity
ini_velocity <- ssa_steady$ini_velocity
steady_thickness <- ssa_steady$current_thickness
ini_thickness <- ssa_steady$ini_thickness
steady_top_surface <- ssa_steady$top_surface
steady_bottom_surface <- ssa_steady$bottom_surface
steady_bed <- ssa_steady$bedrock
steady_friction_coef <- ssa_steady$friction_coef
steady_GL_position <- ssa_steady$grounding_line
x <- ssa_steady$domain
J <- length(x)
# all_u <- ssa_out$all_velocities
# all_H <- ssa_out$all_thicknesses
# all_zs <- ssa_out$all_top_surface

## Plotting

png("./plots/steady_state_matlab/steady_state.png", width = 1000, height = 1000)
par(mfrow = c(2, 2))

# Plot bed and surface elevation -- bed not available for floating ice
plot(domain / 1000, steady_bed,
type = "l", lwd = 2, ylim = c(min(steady_bed), max(steady_top_surface)),
xlab = "Domain (km)", ylab = "Elevation (m)"
) # bed elevation
lines(domain / 1000, steady_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # surface elevation
lines(domain / 1000, steady_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # bottom surface elevation
title("Top and bottom surface elevation")

## Plot GL position
plot(steady_GL_position, 1:length(steady_GL_position), type = "l", xlab = "Domain (km)", ylab = "Iterations")
title("Grounding line position")

## Plot ice thickness
fin_GL_position <- steady_GL_position[length(steady_GL_position)]
plot_region <- 1:J #(fin_GL_position - floor(J / 6)):(fin_GL_position + J / 2) # 500:1200 ##
plot(x[plot_region] / 1000, steady_thickness[plot_region], type = "l", ylab = "Ice thickness (m)")
lines(x[plot_region] / 1000, ini_thickness[plot_region], col = "cyan", lty = 2)
legend("topright",
legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35)
)
title("Ice thickness")

## Plot velocity
plot(domain / 1000, steady_velocity, type = "l", ylim = c(min(steady_velocity), max(steady_velocity)), xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
lines(domain / 1000, ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
legend("topleft",
legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, xpd = TRUE, inset = c(0, -0.35)
)
title("Horizontal velocity")

dev.off()

## Save steady state
if (rewrite_steady_state) {
    ## Save output
    # saveRDS(ssa_out, file = "/home/babv971/SSA_model/EnKF/Output/ssa_steady.rds")
    saveRDS(ssa_out, file = paste("./output/ssa_steady_", data_date, ".rds", sep = ""))
}

## Reduce the ice rigidity to induce dynamics

# ref <- create_ref(domain = ssa_steady$domain,
#       bedrock = ssa_steady$bedrock,
#       friction_coef = ssa_steady$friction_coef,
#       ini_velocity = ssa_steady$current_velocity,
#       ini_thickness = ssa_steady$current_thickness,
#       years = years, steps_per_yr = steps_per_yr,
#       seed = ssa_seed, save_model_output = TRUE,
#       perturb_hardness = TRUE,
#       add_process_noise = add_process_noise_in_ref)

# if (rewrite_reference) {
#     # Save reference data to a file
#     saveRDS(reference, file = paste("./output/reference_", data_date, ".rds", sep = ""))
# }