##### Solve for ice thickness #####

## Notation: u = velocity, H = thickness

solve_thickness <- function(p, velocity, thickness, domain, bed, include_GL = TRUE, fixed_H0 = FALSE,
                            steps_per_yr = 100) {
  
  u <- velocity
  H <- thickness
  
  ## Domain
  x <- domain
  L <- x[length(x)]
  dx <- x[2] - x[1]
  J <- L / dx
  
  ## Convert input velocity into m/s
  # u <- u / secpera
  
  # Surface mass balance rate
  as <- p$as #0.5 / secpera # surface accumulation rate (m/s)
  ab <- p$ab # melt rate (m/s)
  Mg <- as - ab # surface mass balance

  M <- rep(Mg, length(x))
  
  # if (include_GL) {
  #   GL <- find_gl(H, b, rho_i = p$rho_i, rho_w = p$rho_w)
  # }
  M_stag <- 0.5 * (M[2:(J+1)] + M[1:J])
  
  # Time increment
  dt <- p$secpera / steps_per_yr
# NEED TO CHECK IF THE ICE ACTUALLY ACCUMULATES HERE

  # if (transient) {
    # Upwind scheme
    H[2:(J+1)] <- H[2:(J+1)] + dt * M_stag - dt/dx * (u[2:(J+1)] * H[2:(J+1)] - u[1:J] * H[1:J])
  # } else { 
  #   # Steady state
  #   H[2:(J+1)] <- (dx * M_stag + u[1:J] * H[1:J]) / u[2:(J+1)] 
  # }

  # Central difference scheme -- unused
  # H[2:J] <- H[2:J] + dt * M[2:J] - dt/(2 * dx) * (u[3:(J+1)] * H[3:(J+1)] - u[1:(J-1)] * H[1:(J-1)])
  
  # Boundary condition at x = 0
  if (fixed_H0) {
    H[1] <- H0
  } else {
    H[1] <- H[2] #flat thickness at x = 0
  }
  
  H
}