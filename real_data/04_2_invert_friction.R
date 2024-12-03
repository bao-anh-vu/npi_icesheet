## Invert for initial friction
setwd("~/SSA_model/CNN/real_data/")

library(qs)

data_dir <- "./data/"
data_date <- "20241103"


secpera <- 31556926
n <- 3
m <- 1/3
rho <- 910
g <- 9.81

## Domain
## Flowline data
flowline <- readRDS(paste0(data_dir, "/flowline_regrid.rds"))
J <- 2001 # number of grid points
# flowline <- flowline[1:J, ]
flowline_dist <- sqrt((flowline$x[2:J] - flowline$x[1:(J-1)])^2 + (flowline$y[2:J] - flowline$y[1:(J-1)])^2)
flowline_dist <- c(0, cumsum(na.omit(flowline_dist)))

## Velocity data
vel_mat <- qread("./data/velocity/all_velocity_arr.qs")
vel_curr <- vel_mat[, ncol(vel_mat)]
## Smooth the velocity out with loess
vel_curr_smooth <- loess(vel_curr ~ flowline_dist, span = 0.1)
u <- vel_curr_smooth$fitted

## Surface elevation data
surf_elev_mat <- qread("./data/surface_elev/surf_elev_mat.qs")
z <- rowMeans(surf_elev_mat) #, na.rm = T)

## Bedrock data
b <- qread(file = paste0(data_dir, "training_data/bed_sim_steady_state.qs"))

H <- z - b

## Now "invert" for the friction coefficient
x <- flowline_dist
L <- x[length(x)]
  # dx <- x[2] - x[1]
dx <- mean(x[2:length(x)] - x[1:(length(x)-1)]) 
# J <- round(L / dx)

## Ice rigidity
Bg <- 0.4 * 1e6 * 31556926 ^ (1/3)
B <- rep(Bg, length(x))

## Set up some vectors for gamma, beta, W
# gamma <- rep(0, length(x))
# beta <- rep(0, length(x))

## Define a strain rate for regularisation of du/dx
eps_reg <- (1.0 / secpera) / L # strain rate of 1 m/a over length of shelf

# 1. Compute du_dx on a grid and regularise
du_dx <- c( t((u[2] - u[1])/dx), t((u[3:(J+1)] - u[1:(J-1)])/(2*dx)),
            t((u[J+1] - u[J])/dx))

du_dx_reg <- sqrt(du_dx^2 + eps_reg^2) #regularise to avoid division by zero

# 2. Average values of H on a staggered grid (for stability)
H <- c(H, H[J]) ## add a dummy value at the end
H_stag <- 0.5 * (H[2:length(H)] + H[1:(length(H)-1)])
B_stag <- 0.5 * (B[2:(J+1)] + B[1:J])

# and calculate W(x) using the current iteration of du/dx
W <- rep(0, length(x))
W[1:J] <- 2 * B_stag[1:J] * H_stag * abs(du_dx_reg[1:J])^((1/n)-1)
W[J+1] <- W[J] # technically these are W_{J+3/2} and W_{J+1/2}

lhs <- (W[2:(J+1)] * du_dx_reg[2:(J+1)] - W[1:J] * du_dx_reg[1:J])/(dx^2)

# 4. Construct vector b on the rhs
# First compute dz/dx
# z <- c() # top surface elevation
# if (include_GL) {
#     z[1:GL] <- H[1:GL] + b[1:GL]
#     z[(GL+1):(J+1)] <- H[(GL+1):(J+1)] * (1 - rho / rho_w) + z0
# } else {
#     z <- H + b
# }

z <- c(z, z[J]) ## add a dummy value at the end
dz_dx <- c( t((z[2] - z[1])/dx), t((z[3:(J+1)] - z[1:(J-1)])/(2*dx)),
            t((z[J+1] - z[J])/dx))

# then compute the whole rhs

beta <- as.vector(rho * g * H * dz_dx)
rhs <- dx^2 * beta
rhs <- rhs[1:J]

## Basal shear stress
shear_stress <- rhs - lhs
C <- shear_stress / (abs(u)^(m-1) * u)
# C2 <- shear_stress / (u^m)

fric_scale <- 1e6 * secpera^m
C <- C / fric_scale

print("Plotting...")
png(paste0("./plots/friction/friction_inverted.png"), width = 1000, height = 500)
plot(C, type = "l", xlab = "Distance along flowline (km)", ylab = "Friction coefficient (Pa m^-1 s^(1-m))")
# lines(C2, col = "red")
# plot(z, type = "l")
# plot(u, type = "l")
dev.off()
# gamma[1:GL] <- C[1:GL] * abs(u[1:GL])^(m-1)
# gamma[(GL+1):(J+1)] <- 0

