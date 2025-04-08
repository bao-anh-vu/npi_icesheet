setwd("/home/babv971/SSA_model/CNN/alexr/")

library(R.utils)
library(pracma)
library(rootSolve)

source("./source/flowline_eqns.R")
source("./source/create_basal_params.R")

## Flags
use_linear_bed <- F

# Bed topography function
bed <- function(x, params) {
  b <- params$b0 + params$bx * x
  
  xsill <- which(x > params$sill_min & x < params$sill_max)
  xdsill <- which(x >= params$sill_max)
  sill_length <- params$sill_max - params$sill_min
  
  b[xsill] <- params$b0 + (params$bx * params$sill_min) +
              params$sill_slope * (x[xsill] - params$sill_min)
  b[xdsill] <- params$b0 + (params$bx * params$sill_min) +
               params$sill_slope * sill_length +
               params$bx * (x[xdsill] - params$sill_max)
  
  return(b)
}

# Initialize parameters
params <- list()

# Bed parameters
params$b0 <- -100
params$bx <- -1e-3
params$sill_min <- 2000e3
params$sill_max <- 2100e3
params$sill_slope <- 1e-3

# Physical parameters
params$year <- 3600 * 24 * 365
params$Aglen <- 4.227e-25
params$nglen <- 3
params$Bglen <- params$Aglen^(-1 / params$nglen)
params$m <- 1 / params$nglen
params$accum <- 0.28 / params$year
params$C <- 7e6
params$rhoi <- 917
params$rhow <- 1028
params$g <- 9.81

# Scaling parameters
params$hscale <- 1000
params$ascale <- 0.1 / params$year
params$uscale <- (params$rhoi * params$g * params$hscale * params$ascale / params$C)^(1 / (params$m + 1))
params$xscale <- params$uscale * params$hscale / params$ascale
params$tscale <- params$xscale / params$uscale
params$eps <- params$Bglen * ((params$uscale / params$xscale)^(1 / params$nglen)) /
               (2 * params$rhoi * params$g * params$hscale)
params$lambda <- 1 - (params$rhoi / params$rhow)
params$transient <- 0

# Grid parameters
params$tfinal <- 10e3 * params$year
params$Nt <- 1e2
params$dt <- params$tfinal / params$Nt
params$Nx <- 200
params$N1 <- 100
params$sigGZ <- 0.97

sigma1 <- seq(params$sigGZ / (params$N1 + 0.5), params$sigGZ, length.out = params$N1)
sigma2 <- seq(params$sigGZ, 1, length.out = params$Nx - params$N1 + 1)
params$sigma <- c(sigma1, sigma2[-1])

params$sigma_elem <- c(0, 0.5 * (params$sigma[1:(params$Nx - 1)] + params$sigma[2:params$Nx]))
params$dsigma <- diff(params$sigma)

# Override accumulation for steady state
params$accum <- 1 / params$year

# Grounding line position (non-dimensional)
xg <- 200e3 / params$xscale

# # Domain
x <- seq(from = 0, to = xg, length.out = params$Nx)

if (use_linear_bed) {
     hf <- (-bed(xg * params$xscale, params) / params$hscale) / (1 - params$lambda)
     # Bed profile
    b <- -bed(xg * params$sigma * params$xscale, params) / params$hscale
} else {
     # Bed topography
     b <- create_bed(x, variable_bed = F, random_bed = F)

     # Floatation ice thickness at GL (non-dimensional)
     hf <- (-b[which(x == xg)] / params$hscale) / (1 - params$lambda)
}

# Initial ice thickness profile
h <- 1 - (1 - hf) * params$sigma

# Initial velocity profile
u <- 0.3 * (params$sigma_elem)^(1 / 3) + 1e-3

# Convert non-dimensional values to dimensional
h_dim <- h * params$hscale
u_dim <- u * params$uscale * params$year


png("./plots/ini_conds.png")
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Ice thickness plot
plot(params$sigma_elem, h_dim, type = "l", col = "blue", lwd = 2,
     ylim = c(-1000, 1500),
     xlab = "Domain (sigma)", ylab =
      "Ice thickness (m)",
     main = "Initial Conditions")
# lines(params$sigma_elem, b)

# Velocity plot
plot(params$sigma, u_dim, type = "l", col = "red", lwd = 2,
     xlab = "Domain (sigma)", ylab = "Velocity (m/yr)")
dev.off()


## Steady state ##

# Initial guess
huxg0 <- c(h, u, xg)

# Solve
steady_sol <- multiroot(
  f = function(huxg) flowline_eqns(huxg, b, params),
  start = huxg0,
  maxiter = 1000,
  rtol = 1e-8
)


# Unpack solution
h <- steady_sol$root[1:params$Nx]
u <- steady_sol$root[(params$Nx + 1):(2 * params$Nx)]
xg <- steady_sol$root[2 * params$Nx + 1]

cat(sprintf("Steady state grounding line: %.2f km\n", xg * params$xscale / 1e3))

# Create dimensional spatial grid
x_dim <- params$sigma * xg * params$xscale / 1e3  # in km

# Plot ice thickness
# Create dimensional spatial grid for sigma elements (ice thickness h)
x_elem <- params$sigma_elem * xg * params$xscale / 1e3  # in km

# Bedrock topography at same points
bed_profile <- bed(x_elem * 1e3, params)  # in meters

# Plot h and bedrock b together

png("./plots/steady_state.png", width = 1000, height = 800)

par(mfrow = c(2, 1))

plot(x_elem, h * params$hscale,
     type = "l", lwd = 2, col = "steelblue",
     ylim = range(c(bed_profile, h * params$hscale)),
     xlab = "Distance from ice divide (km)",
     ylab = "Elevation (m)",
     main = "Steady-State Ice Thickness and Bedrock")
lines(x_elem, bed_profile, lwd = 2, col = "sienna")
legend("topright", legend = c("Ice thickness (h)", "Bedrock (b)"),
       col = c("steelblue", "sienna"), lwd = 2)

# Plot ice velocity
plot(x_dim, u * params$uscale,
     type = "l", lwd = 2, col = "firebrick",
     xlab = "Distance from ice divide (km)",
     ylab = "Ice velocity u(x) (m/yr)",
     main = "Steady-State Ice Velocity")

abline(v = xg * params$xscale / 1e3, col = "gray40", lty = 2)  # grounding line marker

dev.off()

