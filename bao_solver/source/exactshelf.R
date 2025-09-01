exactshelf <- function(x, L = 200e3, M0 = NULL, Hg = 500, ug = NULL) {
  # Physical parameters from MISMIP
  param <- list(
    secpera = 31556926,
    n = 3.0,
    rho = 900.0,
    rhow = 1000.0,
    g = 9.8,
    A = 1.4579e-25
  )

  # Set defaults if needed
  if (is.null(M0)) M0 <- 0.3 / param$secpera
  if (is.null(ug)) ug <- 50 / param$secpera

  # Derived constants
  n <- param$n
  r <- param$rho / param$rhow
  Cs <- param$A * (0.25 * param$rho * param$g * (1 - r))^n
  qg <- ug * Hg

  # Allocate and compute only where x in [0, L]
  u <- rep(0, length(x))
  H <- rep(0, length(x))
  inside <- which(x >= 0 & x <= L)

  u[inside] <- (ug^(n + 1) + (Cs / M0) *
                  ((M0 * x[inside] + qg)^(n + 1) - qg^(n + 1)))^(1 / (n + 1))
  H[inside] <- (M0 * x[inside] + qg) / u[inside]

  return(list(uexact = u, H = H))
}
