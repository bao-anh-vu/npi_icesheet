setwd("/home/babv971/SSA_model/CNN/bao_solver/")

rm(list = ls())


# Parameters
L <- 1000000
nx <- 300
dx <- L / (nx - 1)
x <- seq(0, L, length.out = nx)

# Bed profile: linear decline below sea level
b <- -500 + 0.1 * x / 1000  # bed drops from -500 m to -400 m
rho <- 917
rho_w <- 1028
g <- 9.81
A <- 1e-16
n <- 3
a <- rep(0.3 / (365 * 24 * 3600), nx)  # accumulation in m/s

H <- rep(100, nx)
u <- rep(0, nx)

dt <- 0.01 * 365 * 24 * 3600
nt <- 1000

# Function: check flotation condition
is_floating <- function(H, b, rho, rho_w) {
  return(H < -b * rho_w / rho)
}

# SSA with basal drag: zero if floating
solve_velocity <- function(H, b, dx, A, n, rho, g, rho_w) {
  u <- rep(0, length(H))
  for (i in 2:(length(H) - 1)) {
    dHdx <- (H[i+1] - H[i-1]) / (2 * dx)
    h <- H + b
    dhdx <- (h[i+1] - h[i-1]) / (2 * dx)
    eta <- 0.5 * A^(-1/n) * H[i] * abs(dHdx)^(1/n - 1)
    
    # Basal drag if grounded, zero if floating
    floating <- is_floating(H[i], b[i], rho, rho_w)
    tau_b <- if (floating) 0 else 100 * u[i]  # simple linear friction law

    u[i] <- (rho * g * H[i] * dhdx - tau_b) / (2 * eta)
  }
  u[1] <- 0
  u[length(u)] <- u[length(u) - 1]
  return(u)
}

# Time loop
for (t in 1:nt) {
  u <- solve_velocity(H, b, dx, A, n, rho, g, rho_w)
  q <- u * H
  dqdx <- c(0, diff(q) / dx)
  H_new <- H + dt * (a - dqdx)
  H_new[H_new < 0] <- 0
  H <- H_new
  
  if (t %% 100 == 0) {
    floating_mask <- is_floating(H, b, rho, rho_w)
    grounding_line <- max(x[!floating_mask], na.rm = TRUE)
    

    png(paste0("output_", t, ".png"), width = 800, height = 600)
    plot(x / 1000, H + b, type = 'l', col = "blue", ylim = c(-600, 1500),
         xlab = "Distance (km)", ylab = "Surface elevation (m)",
         main = paste("Time =", round(t * dt / (365 * 24 * 3600)), "years"))
    lines(x / 1000, b, col = "black")
    abline(v = grounding_line / 1000, col = "red", lty = 2)
    legend("topright", legend = c("Surface", "Bed", "Grounding line"),
           col = c("blue", "black", "red"), lty = c(1,1,2))
    dev.off()
  }
}
