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

## Parameters
secpera <- 31556926 # seconds per annum
n <- 3.0 # exponent in Glen's flow law
m <- 1 / n # friction exponent

## Ice shelf specifications
L <- 800e3 # length of domain

## Boundary conditions at the ice divide (x = 0)
x0 <- 0
u0 <- 100 / secpera # (m/s)
H0 <- 2000 # (m)

# Discretise domain
J <- 2000 # number of grid points
dx <- L / J # grid spacing
x <- seq(x0, L, dx) # grid point locations

## Params

### Friction
C <- create_fric_coef(x, L) * 1e6 * (secpera)^m

### Bed
b <- create_bed(x)

ssa_steady <- create_steady_state(use_stored_steady_state = F,
                        add_process_noise_in_ref = F, 
                        rewrite_steady_state = F,
                        domain = x, 
                        friction_coef = C,
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
par(mfrow = c(2, 2))

# Plot bed and surface elevation -- bed not available for floating ice
plot(x / 1000, steady_bed,
type = "l", lwd = 2, xlim = c(0, 800), ylim = c(min(steady_bed), max(steady_top_surface)),
xlab = "Domain (km)", ylab = "Elevation (m)"
) # bed elevation
lines(x / 1000, steady_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # surface elevation
lines(x / 1000, steady_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # bottom surface elevation
title("Top and bottom surface elevation")

## Plot GL position
plot(steady_GL_position, 1:length(steady_GL_position), type = "l", xlab = "Domain (km)", ylab = "Iterations")
title("Grounding line position")

## Plot ice thickness
fin_GL_position <- steady_GL_position[length(steady_GL_position)]
plot_region <- (fin_GL_position - floor(J / 6)):(fin_GL_position + J / 2) # 500:1200 ##
plot(x[plot_region] / 1000, steady_thickness[plot_region], type = "l", ylab = "Ice thickness (m)")
lines(x[plot_region] / 1000, ini_thickness[plot_region], col = "cyan", lty = 2)
legend("topright",
legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35)
)
title("Ice thickness")

## Plot velocity
plot(x / 1000, steady_velocity, type = "l", ylim = c(min(steady_velocity), max(steady_velocity)), xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
lines(x / 1000, ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
legend("topleft",
legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, xpd = TRUE, inset = c(0, -0.35)
)
title("Horizontal velocity")

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