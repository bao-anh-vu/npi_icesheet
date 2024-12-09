## Function to propagate
propagate <- function(state, domain, steps_per_yr, transformation = "log10") {
  
  thickness <- state[1:J]
  bed <- state[(J+1):(2*J)] 
  alpha <- state[(2*J+1):(3*J)]
  velocity <- state[(3*J+1):(4*J)]
  
  friction <- NULL
  if (transformation == "log") {
    friction <- exp(alpha)
  } else if (transformation == "sqrt") {
    friction <- alpha^2
  } else {
    friction <- 10^alpha
  }
  
  prev_velocity <- c()
  
  for (i in 1:steps_per_yr) {
    prev_velocity <- velocity
    
    thickness <- solve_thickness(velocity = prev_velocity,
                                 thickness = thickness, domain = domain,
                                 bed = bed, steps_per_yr = steps_per_yr)
    
    browser()
    ## need to somehow save both the previous and the current velocity
    velocity <- as.vector(solve_velocity(prev_velocity = prev_velocity, 
                                         thickness = thickness, # should use the mean of the ice thickness as well?
                                         domain = domain,
                                         bed = bed,
                                         friction = friction,
                                         increase_hardness = F))
    
  }
  
  out <- c(thickness, bed, alpha, prev_velocity)
  # out <- c(thickness, bed, friction, velocity)
  
  return(out)
}