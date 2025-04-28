##### Solve for ice thickness #####

## Notation: u = velocity, H = thickness

solve_thickness <- function(velocity, thickness, domain, bed, include_GL = TRUE, fixed_H0 = FALSE,
                            M_variable = FALSE, M_bueler = FALSE,
                            secpera = 31556926, steps_per_yr = 100) {
  
  u <- velocity
  H <- thickness
  
  ## Domain
  x <- domain
  L <- x[length(x)]
  dx <- x[2] - x[1]
  J <- L / dx
  
  ## Convert input velocity into m/s
  u <- u / secpera
  
  # Surface mass balance rate
  if (M_bueler) {
    Mg <- -4.290 / secpera #(m/s) # SMB at GL
  } else {
    as <- 0.5 / secpera # surface accumulation rate (m/s)
    ab <- 0 # melt rate (m/a)
    Mg <- as - ab # surface mass balance
  }
  M <- rep(Mg, length(x))
  
  if (include_GL) {
    GL <- find_gl(H, bed)
    
    if (M_variable) {
      M[1:GL] <- a * (H[1:GL] - Hela)
    }
  }
  M_stag <- 0.5 * (M[2:(J+1)] + M[1:J])
  
  # Time increment
  dt <- secpera / steps_per_yr
  
  # Upwind scheme
  H[2:(J+1)] <- H[2:(J+1)] + dt * M_stag - dt/dx * (u[2:(J+1)] * H[2:(J+1)] - u[1:J] * H[1:J])
  
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