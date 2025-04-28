# testshelf <- function(J = 300, L = 200e3) {
setwd("~/SSA_model/mccarthy/rfiles")

rm(list = ls())

library(matrixStats)

source("ssaflowline_bao.R") # function to solve for the velocity
source("flowline_bao.R")
# source("ssainit.R")
source("find_gl.R")
source("create_params.R") # function to create the bedrock profile

J <- 2000
L <- 800e3

# Physical parameters
param <- list(
    secpera = 31556926,
    n = 3.0,
    rho = 900.0,
    rhow = 1000.0,
    g = 9.8,
    A = 1.4579e-25
)

param$m <- 1 / param$n
param$L <- L
dx <- param$L / J
x <- seq(0, param$L, length.out = J + 1)

# Default ice stream-shelf properties
# Changed these boundary conditions from grounding line conditions to ice divide instead
# Hg <- 500 # thickness at the grounding line
# ug <- 50 / param$secpera # velocity at the grounding line
# u_bc <- 0 / param$secpera # velocity at the ice divide (x1 = 0)

H <- 2100 - 2000 / L * x
# b <- ifelse(x <= 450e3, -1100 + x / 1e3, -650 - 5 * (x / 1e3 - 450))
b <- create_bed(x)

GL <- find_gl(H, b, rho_i = param$rho, rho_w = param$rhow)

# Compute bottom surface elevation for the shelf
r <- param$rho / param$rhow
zb_shelf <- -r * H[(GL+1):(J+1)]

zb <- c(b[1:GL], zb_shelf)
zs <- H + zb # surface elevation (z_s in my notation)
# zs[(GL+1):(J+1)] <- (1 - param$rho/param$rhow) * H[(GL+1):(J+1)] # for floating ice

# Basal drag
C <- rep(0, length(x)) # for floating ice
C[1:GL] <- 7e6 * (param$secpera^param$m) # for grounded ice; maybe try something linear

# Initial velocity (modify later to use ssainit())
# ug <- 50 / param$secpera # velocity at the grounding line
# gamma <- (0.25 * param$A^(1 / param$n) * (1 - param$rho / param$rhow) * param$rho * param$g * H[J + 1])^param$n
# u0_shelf <- ssainit(param, x_shelf, ug, drivestress = 0, gamma, initchoice = 1) # drivestress is not used here so set to 0
# u0_grounded <- seq(0, ug, length.out = J + 1 - length(x_shelf))
# u0 <- c(u0_grounded, u0_shelf)

u0 <- 0.0001 * x / param$secpera
u_bc <- 0 # boundary condition at the ice divide (x1 = 0)
# u0 <- ssainit(p, x, beta, gamma, initchoice)

## Plot initial conditions
png("ini_conds.png", width = 800, height = 800)
par(mfrow = c(2, 1), mar = c(5, 5, 2, 1))

plot(x / 1000, zs,
    type = "l", col = "blue", lwd = 3,
    xlab = "x (km)", ylab = "thickness (m)", ylim = range(c(H, b))
)
lines(x / 1000, zb, col = "blue", lwd = 3)
lines(x / 1000, b, col = "red", lwd = 3)
legend("topleft",
    legend = c("Top surface", "Bedrock"),
    col = c("blue", "red"), lty = c(1, 1, NA),
    bty = "n", lwd = 3
)

plot(x / 1000, u0 * param$secpera,
    type = "l", col = "blue", lwd = 3,
    xlab = "x (km)", ylab = "velocity (m/a)", ylim = range(u0) * param$secpera
)
legend("topleft",
    legend = c("Initial guess", "Numerical"),
    col = c("blue", "red"), lty = c(1, 1, NA), pch = c(NA, NA, 1),
    bty = "n", lwd = 3
)
#   }
dev.off()

# Numerical solution
result <- solve_velocity(param, J, H, b, C, u0, u_bc)
unum <- result$u
u0 <- result$u0
zs <- result$zs
H <- result$H
zb <- zs - H

zb[(GL + 1):(J + 1)] <- zs[(GL + 1):(J + 1)] - H[(GL + 1):(J + 1)] # for floating ice

# Compute errors
# averr <- sum(abs(unum - uexact)) / (J + 1)
# maxerr <- max(abs(unum - uexact))

# Plot if not called for return
#   if (interactive()) {
png("steady_state.png", width = 800, height = 800)
par(mfrow = c(2, 1), mar = c(5, 5, 2, 1))

plot(x / 1000, zs,
    type = "l", col = "blue", lwd = 3,
    xlab = "x (km)", ylab = "thickness (m)", ylim = range(c(H, b))
)
lines(x / 1000, zb, col = "blue", lwd = 3)
lines(x / 1000, b, col = "red", lwd = 3)
legend("topleft",
    legend = c("Ice surface", "Bedrock"),
    col = c("blue", "red"), lty = c(1, 1, NA),
    bty = "n", lwd = 3
)

plot(x / 1000, u0 * param$secpera,
    type = "l", col = "blue", lwd = 3,
    xlab = "x (km)", ylab = "velocity (m/a)", ylim = range(unum) * param$secpera
)
# lines(x / 1000, uexact * param$secpera, col = "red", lwd = 3)
points(x / 1000, unum * param$secpera, lwd = 3, type = "l", col = "red")
legend("topleft",
    legend = c("Initial guess", "Numerical"),
    col = c("blue", "red"), lty = c(1, 1, NA), pch = c(NA, NA, 1),
    bty = "n", lwd = 3
)
#   }
dev.off()

# cat(sprintf("Average error = %.4e m/s\n", averr))
# cat(sprintf("Maximum error = %.4e m/s\n", maxerr))

#   return(list(averr = averr, maxerr = maxerr))
# }
