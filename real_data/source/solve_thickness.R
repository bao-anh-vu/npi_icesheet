##### Solve for ice thickness #####

## Notation: u = velocity, H = thickness

solve_thickness <- function(velocity, thickness, domain, bed, include_GL = TRUE, fixed_H0 = FALSE,
                            M_variable = FALSE, M_bueler = FALSE,
                            secpera = 31556926, steps_per_yr = 100,
                            # use_relaxation = FALSE, #relax_thickness = NULL,
                            relax_rate = 0, as = 0.5, ab = 0) {
  
  u <- velocity
  H <- thickness
  
  ## Domain
  x <- domain
  L <- x[length(x)]
  dx <- mean(x[2:length(x)] - x[1:(length(x)-1)]) #x[2] - x[1]
  J <- round(L / dx)
  
  ## Convert input velocity into m/s
  u <- u / secpera
  
  # Surface mass balance rate
  # if (M_bueler) {
  #   Mg <- -4.290 / secpera #(m/s) # SMB at GL
  # } else {
  
# browser()
# plot(x / 1000, as, type = "l", col = "blue", ylim = c(-2.5, 1.2),
#       xlab = "Domain (km)", ylab = "Accumulation rate (m/a)")
# lines(x / 1000, relax_rate, col = "red")
# lines(x / 1000, as + relax_rate, col = "purple")
# abline(v = domain[GL] / 1000, lty = 2, lwd = 2)


    # as <- as / secpera # surface accumulation rate (m/s)
    # ab <- ab / secpera # melt rate (m/s)
    Mg <- (as - ab) / secpera# surface mass balance
  # }

  if (length(Mg) == 1) {
    M <- rep(Mg, length(x))
  } else {
    M <- Mg
  }
  
  if (include_GL) {
    GL <- gl_migrate(H, bed)
    
    if (M_variable) {
      M[1:GL] <- a * (H[1:GL] - Hela)
    }
  }
  M_stag <- 0.5 * (M[2:(J+1)] + M[1:J]) # compute on staggered grid

# if (use_relaxation) {
  relax_rate <- relax_rate / secpera

  if (length(relax_rate) == 1) {
    relax_rate <- rep(relax_rate, J+1)
  }

  relax_rate_stag <- 0.5 * (relax_rate[2:(J+1)] + relax_rate[1:J]) # compute on staggered grid
  # relax_term <- 1e-3 / secpera * (relax_thickness - H) # should probably stagger this too?
  # relax_term <- relax_rate
  # M_stag <- M_stag + (u[2:(J+1)] * H[2:(J+1)] - u[1:J] * H[1:J])/dx + relax_term[1:J]
  M_stag <- M_stag + relax_rate_stag
# }

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