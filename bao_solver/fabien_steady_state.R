### Main file ###

setwd("/home/babv971/SSA_model/CNN/ssa_solver/")

rm(list = ls())

# library("parallel")
library("Matrix")
library("qlcMatrix")
library("fastmatrix")
library("expm")
library("R.utils")
# library("sp")
library("fields")
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
source("./source/solve_ssa_nl.R")
source("./source/solve_velocity_azm.R")
source("./source/solve_thickness.R")

## Presets
data_date <- "20250408" # "20230518"
rewrite_steady_state <- F

## Parameters
secpera <- 31556926 # seconds per annum
n <- 3.0 # exponent in Glen's flow law
m <- 1 / n # friction exponent

## Ice shelf specifications
L <- 800e3 # length of domain

# Discretise domain
J <- 2000 # number of grid spaces
dx <- L / J # grid spacing
x <- seq(0, L, dx) # grid point locations

## Params
params <- list(
  secpera = 31556926, # seconds per annum
  n = 3.0, # exponent in Glen's flow law
  rho_i = 917.0, # ice density
  rho_w = 1028.0, # sea water density
  # r <- rho / rho_w # define this ratio for convenience
  g = 9.81, # gravity constant
  # A = 4.227e-25, #1.4579e-25, # flow rate parameter
  as = 0.5, # surface accumulation rate,
  ab = 0 # basal melt rate,
)

params$m <- 1 / params$n
params$B <- 0.4 * 1e6 * params$secpera^params$m

### Friction
C <- create_fric_coef(x, L) * 1e6 * (secpera)^m

### Bed
b <- create_bed(x)

## Boundary conditions at the ice divide (x = 0)
u0 <- 0 # 100 / secpera # (m/s)
H0 <- 2000 # (m)

### Initial thickness
h <- H0 - (H0 - 0) / (L - 0) * x

### Initial velocity
u <- u0 + 0.001 / params$secpera * x

# ssa_steady <- create_steady_state(use_stored_steady_state = F,
#                         add_process_noise_in_ref = F,
#                         rewrite_steady_state = F,
#                         phys_params = params,
#                         domain = x,
#                         ini_thickness = h,
#                         ini_velocity = u,
#                         friction_coef = C,
#                         bedrock = b)

ssa_steady <- solve_ssa_nl(
  domain = x,
  ini_thickness = h,
  ini_velocity = u,
  bedrock = b,
  friction_coef = C,
  tol = 1e-02, # m/yr
  years = 2000, # default is 4000 years but could increase further
  steps_per_yr = 52,
  phys_params = params,
  perturb_hardness = FALSE,
  add_process_noise = FALSE
)

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
# x <- ssa_steady$domain
# J <- length(x)

# all_u <- ssa_out$all_velocities
# all_H <- ssa_out$all_thicknesses
# all_zs <- ssa_out$all_top_surface

## Plotting

png("./plots/steady_state/steady_state.png", width = 1000, height = 800)
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
plot_region <- (fin_GL_position - floor(J / 10)):(fin_GL_position + J / 2) # 500:1200 ##
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

dev.off()

## Save steady state
if (rewrite_steady_state) {
  saveRDS(ssa_out, file = paste("./output/ssa_steady_", data_date, ".rds", sep = ""))
}

## Reduce the ice rigidity to induce dynamics, and run for another 20 years post-steady state
params$B <- 0.25 * 1e6 * params$secpera^params$m # reduce from 0.4 to 0.25

### For this part I also add process noise

### Process noise parameters (for ice thickness)
exp_cov <- function(d, l) {
  return(exp(-3 * d / l))
}

ones <- rep(1, length(x))
D <- rdist(x)
l <- 50e3
R <- exp_cov(D, l)


# R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
L <- t(chol(R))
L <- as(L, "dgCMatrix")
process_noise_info <- list(corrmat_chol = L, length_scale = l)

post_ss <- solve_ssa_nl(
  domain = ssa_steady$domain,
  bedrock = ssa_steady$bedrock,
  friction_coef = ssa_steady$friction_coef,
  phys_params = params,
  ini_velocity = ssa_steady$current_velocity,
  ini_thickness = ssa_steady$current_thickness,
  years = 20, steps_per_yr = 52,
  add_process_noise = T,
  process_noise_info = process_noise_info
)

## Extract individual components
post_steady_velocity <- post_ss$current_velocity
post_steady_ini_velocity <- post_ss$ini_velocity
post_steady_thickness <- post_ss$current_thickness
post_steady_ini_thickness <- post_ss$ini_thickness
post_steady_top_surface <- post_ss$top_surface
post_steady_bottom_surface <- post_ss$bottom_surface
post_steady_bed <- post_ss$bedrock
post_steady_friction_coef <- post_ss$friction_coef
post_steady_GL_position <- post_ss$grounding_line
# x <- post_ss$domain
# J <- length(x)

## Plotting

png("./plots/steady_state/post_steady_20yrs.png", width = 1000, height = 800)
par(mfrow = c(2, 2))

# Plot bed and surface elevation -- bed not available for floating ice
plot(x / 1000, post_steady_bed,
  type = "l", lwd = 2, xlim = c(0, 800), ylim = c(min(post_steady_bed), max(post_steady_top_surface)),
  xlab = "Domain (km)", ylab = "Elevation (m)"
) # bed elevation
lines(x / 1000, post_steady_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # surface elevation
lines(x / 1000, post_steady_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # bottom surface elevation
title("Top and bottom surface elevation")

## Plot GL position
plot(post_steady_GL_position, 1:length(post_steady_GL_position), type = "l", xlab = "Domain (km)", ylab = "Iterations")
title("Grounding line position")

## Plot ice thickness
fin_GL_position <- post_steady_GL_position[length(post_steady_GL_position)]
plot_region <- (fin_GL_position - floor(J / 10)):(fin_GL_position + J / 2) # 500:1200 ##
plot(x[plot_region] / 1000, post_steady_thickness[plot_region], type = "l", ylab = "Ice thickness (m)")
lines(x[plot_region] / 1000, post_steady_ini_thickness[plot_region], col = "cyan", lty = 2)
legend("topright",
  legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
  lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35)
)
title("Ice thickness")

## Plot velocity
plot(x / 1000, post_steady_velocity, type = "l", ylim = c(min(post_steady_velocity), max(post_steady_velocity)), xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
lines(x / 1000, post_steady_ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
legend("topleft",
  legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
  lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, xpd = TRUE, inset = c(0, -0.35)
)
title("Horizontal velocity")

dev.off()

