initialise_ice_thickness <- function(domain, surface_obs, bed, 
                                     n_sims,
                                     rho = 910, rho_w = 1028, 
                                     condsim_shelf = FALSE) {

  # Ne <- n_sims
  # 
  # if (is.null(ncol(bed))) {
  #   n_beds <- 1
  # } else {
  #   n_beds <- ncol(bed)
  # }

  # Calculate ice thickness from z and b
  simulated_thickness <- matrix(0, nrow = J, ncol = n_sims)
  thickness_ens <- list()
  
  for (i in 1:n_sims) {
    h <- surface_obs - bed #s[, j]
    
    ## Calculate GL position for each ensemble member
    GL <- gl_migrate2(surface_obs, bed) #s[, i])
    # cat("GL: ", GL, "\n")
    
    ## Recalculate ice thickness from GL onwards (ice shelf)
    shelf_region <- GL:length(h)
    h[shelf_region] <- 1 / (1 - rho / rho_w) * surface_obs[shelf_region]
  
    simulated_thickness[, i] <- h
  }

  return(simulated_thickness)
  
}

# initialise_ice_thickness <- function(domain, surface_obs, bed, 
#                                      n_sims,
#                                      rho = 910, rho_w = 1028, 
#                                      condsim_shelf = FALSE) {

#   # Ne <- n_sims
#   # 
#   # if (is.null(ncol(bed))) {
#   #   n_beds <- 1
#   # } else {
#   #   n_beds <- ncol(bed)
#   # }

#   # Calculate ice thickness from z and b
#   simulated_thickness <- matrix(0, nrow = J, ncol = n_sims)
#   thickness_ens <- list()
  
#   for (i in 1:n_sims) {
#     h <- surface_obs - bed #s[, j]
    
#     ## Calculate GL position for each ensemble member
#     GL <- gl_migrate2(surface_obs, bed) #s[, i])
#     # cat("GL: ", GL, "\n")
    
#     ## Recalculate ice thickness from GL onwards (ice shelf)
#     shelf_region <- GL:length(h)
#     h[shelf_region] <- 1 / (1 - rho / rho_w) * surface_obs[shelf_region]
    
  
#     if (condsim_shelf) { # use conditional simulation on the ice shelf only
#       df <- data.frame(x = domain[shelf_region], h = h[shelf_region])
#       h.fit <- loess(h ~ x, data = df, span = 0.2)
      
#       # h.sd <- seq(10, 1, length.out = length(domain))
#       # shelf.sd <- h.sd[shelf_region]
#       shelf.sd <- pmin(0.02 * h.fit$fitted, 20)
      
      
#       simulated_shelf <- cond_sim(mu = h.fit$fitted, sd = shelf.sd, #pmin(0.02 * h.fit$fitted, 10),
#                                   taxis = domain[shelf_region], l = 50e3, nsims = 1)
#       simulated_thickness[, i] <- c(h[-shelf_region], simulated_shelf[, 1])
      
#     } else { # use conditional simulation on the whole domain
#       ## Need to do this separately for grounded vs floating ice
#       df_float <- data.frame(x = domain[shelf_region], h = h[shelf_region])
#       h.fit_float <- loess(h ~ x, data = df_float, span = 0.2)
      
#       df_ground <- data.frame(x = domain[-shelf_region], h = h[-shelf_region])
#       h.fit_ground <- loess(h ~ x, data = df_ground, span = 0.05)
      
#       h.sd <- c(rep(20, GL), rep(5, length(h) - GL))
#       # h.sd <- pmin(0.02 * c(h.fit_float$fitted, h.fit_ground$fitted), 50)
      
#       sim_float <- cond_sim(mu = h.fit_float$fitted, sd = h.sd[shelf_region],
#                             taxis = domain[shelf_region], l = 50e3, nsims = 1)
      
      
#       sim_ground <- cond_sim(mu = h.fit_ground$fitted, sd = h.sd[-shelf_region],
#                              taxis = domain[-shelf_region], l = 5e3, nsims = 1)
      
#       simulated_thickness[, i] <- c(sim_ground[, 1], sim_float[, 1])#simulated_h[, 1]
#     }
#   }

#   return(simulated_thickness)
  
# }