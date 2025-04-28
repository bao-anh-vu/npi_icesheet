setwd("~/SSA_model/mccarthy/rfiles")

source("flowline.R")
# testflowline <- function(test = 2, J = 300) {
  
  J <- 300
  test <- 2

  L <- 1
  dx <- L / J
  x <- seq(0, L, length.out = J + 1)
  xstag <- seq(dx/2, L + dx/2, by = dx)

  if (test == 1) {
    W <- rep(1, length(xstag))
    alpha <- rep(0, length(x))
    beta <- rep(0, length(x))
    gamma <- 1 / L
    uexact <- x / L

  } else if (test == 2) {
    W <- rep(2, length(xstag))
    alpha <- rep(1, length(x))
    uexact <- sin(x)
    gamma <- cos(1)
    beta <- -3 * sin(x)

  } else if (test == 3) {
    W <- 2 - xstag
    alpha <- rep(1, length(x))
    uexact <- sin(x)
    gamma <- cos(1)
    beta <- -cos(x) - (3 - x) * sin(x)

  } else {
    stop("Only test == 1, 2, or 3 are allowed.")
  }

  # Run flowline
  unum <- flowline(L, J, gamma, W, alpha, beta, ug = 0)
  maxerr <- max(abs(unum - uexact))

  # Plot if not being called inside another function
  plot(x, uexact, type = "l", col = "blue", lwd = 2,
       ylab = "u", xlab = "x", main = paste("Test", test, "with J =", J))
  points(x, unum, col = "red", pch = 1)
  legend("topright", legend = c("Exact", "Numerical"), col = c("blue", "red"),
         lty = c(1, NA), pch = c(NA, 1), bty = "n")

  cat(sprintf("Maximum absolute error = %.4e\n", maxerr))
  return(maxerr)
# }
