## Observation operator ##
obs_operator <- function(state, domain, transformation = "log10", rho = 910, rho_w = 1028) {
  h <- state[1:J] # ice thickness
  b <- state[(J+1):(2*J)] # bed
  alpha <- state[(2*J + 1):(3*J)] # log10(friction)
  velocity <- state[(3*J + 1):(4*J)] # velocity
  
  GL <- gl_migrate(h, b)
  
  ## Calculate surface obs
  z <- rep(0, J)
  z[1:GL] <- h[1:GL] + b[1:GL]
  z[(GL+1):J] <- (1 - rho/rho_w) * h[(GL+1):J]
  
  ## Calculate surface velocity
  C <- NULL
  if (transformation == "log") {
    C <- exp(alpha)
  } else if (transformation == "sqrt") {
    C <- alpha^2
  } else {
    C <- 10^(alpha)
  }
  
  u <- as.vector(solve_velocity(prev_velocity = velocity, h, domain = domain, 
                                bed = b, friction = C,
                                perturb_hardness = TRUE))
  
  obs <- c(z, u)
  
  return(obs)
}