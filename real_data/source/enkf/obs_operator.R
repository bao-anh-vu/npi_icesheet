## Observation operator ##
obs_operator <- function(state, domain, #missing_pattern = NULL,
                        phys_params,
                        transformation = "log") {
  h <- state[1:J] # ice thickness
  b <- state[(J+1):(2*J)] # bed
  alpha <- state[(2*J + 1):(3*J)] # log10(friction)
  velocity <- state[(3*J + 1):(4*J)] # velocity
  
  GL <- gl_migrate(h, b)
  
  ## Calculate surface obs
  z <- rep(0, J)
  z[1:GL] <- h[1:GL] + b[1:GL]
  z[(GL+1):J] <- (1 - phys_params$rho_i/phys_params$rho_w) * h[(GL+1):J]

  ## Calculate surface velocity
  C <- NULL
  if (transformation == "log") {
    C <- exp(alpha)
  } else if (transformation == "sqrt") {
    C <- alpha^2
  } else {
    C <- 10^(alpha)
  }
  
  fric_scale <- 1e6 * phys_params$secpera^(1 / phys_params$n)
  C <- C * fric_scale

  u <- as.vector(solve_velocity(prev_velocity = velocity, 
                                thickness = h, 
                                domain = domain, 
                                bed = b, friction = C,
                                phys_params = phys_params))

  obs <- c(z, u)

  # if (!is.null(missing_pattern)) {
  #   obs[missing_pattern] <- 0 
  # }

  return(obs)
}