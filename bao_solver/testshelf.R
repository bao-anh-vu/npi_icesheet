# testshelf <- function(J = 300, L = 200e3) {
setwd("~/SSA_model/mccarthy/rfiles")

rm(list = ls())

source("exactshelf.R")
source("ssaflowline.R")
source("flowline.R")
source("ssainit.R")

J <- 300
L <- 200e3

  # Physical parameters
  param <- list(
    secpera = 31556926,
    n = 3.0,
    rho = 900.0,
    rhow = 1000.0,
    g = 9.8,
    A = 1.4579e-25
  )

  param$L <- L
  dx <- param$L / J
  x <- seq(0, param$L, length.out = J + 1)

  # Default shelf properties
  Hg <- 500
  ug <- 50 / param$secpera
  M0 <- 0.3 / param$secpera

  # Get exact solution
  shelf <- exactshelf(x, param$L, M0, Hg, ug)  # returns list(uexact, H)
  uexact <- shelf$uexact
  H <- shelf$H

  # Compute flotation bed
  r <- param$rho / param$rhow
  b <- -r * H

  # No basal drag for shelves
  param$C <- 0

  # Numerical solution
  result <- ssaflowline(param, J, H, b, ug, initchoice = 1)
  unum <- result$u
  u0 <- result$u0
  zs <- result$zs


  # Compute errors
  averr <- sum(abs(unum - uexact)) / (J + 1)
  maxerr <- max(abs(unum - uexact))

  # Plot if not called for return
  png("testshelf.png", width = 800, height = 800)
    par(mfrow = c(2, 1), mar = c(5, 5, 2, 1))

    plot(x / 1000, zs, type = "l", col = "blue", lwd = 3,
         xlab = "x (km)", ylab = "thickness (m)", ylim = range(c(H, b)))
    lines(x / 1000, b, col = "red", lwd = 3)
    legend("topleft", legend = c("Top surface", "Bedrock"),
           col = c("blue", "red"), lty = c(1, 1, NA), 
           bty = "n", lwd = 3)

    plot(x / 1000, u0 * param$secpera, type = "l", col = "blue", lwd = 3,
         xlab = "x (km)", ylab = "velocity (m/a)", ylim = range(c(uexact, unum)) * param$secpera)
    lines(x / 1000, uexact * param$secpera, col = "red", lwd = 3)
    points(x / 1000, unum * param$secpera, col = "green", pch = 1, cex = 1.2)
    legend("topleft", legend = c("Initial guess", "Exact solution", "Numerical"),
           col = c("blue", "red", "green"), lty = c(1, 1, NA), pch = c(NA, NA, 1),
           bty = "n", lwd = 3)
  dev.off()
  
  cat(sprintf("Average error = %.4e m/s\n", averr))
  cat(sprintf("Maximum error = %.4e m/s\n", maxerr))

#   return(list(averr = averr, maxerr = maxerr))
# }
