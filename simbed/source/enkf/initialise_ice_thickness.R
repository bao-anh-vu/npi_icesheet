initialise_ice_thickness <- function(domain, surface_obs, bed, 
                                     n_sims,
                                     rho = 910, rho_w = 1028, 
                                     condsim_shelf = FALSE) {

  # Calculate ice thickness from z and b
  simulated_thickness <- matrix(0, nrow = J, ncol = n_sims)
  thickness_ens <- list()
  
  ## Subtract bed from observed surface elevation to get ice thickness
  h <- surface_obs - bed #s[, j]

  ## Recalculate ice thickness from GL onwards (ice shelf)
  ## Since the thickness of floating ice needs to be determined 
  ## based on the flotation condition
  GL <- gl_migrate2(surface_obs, bed) # GL position
  shelf_region <- GL:length(h)
  h[shelf_region] <- 1 / (1 - rho / rho_w) * surface_obs[shelf_region]
  # h_shelf <- 1 / (1 - rho / rho_w) * surface_obs[shelf_region]
  # h <- c(h[1:(GL - 1)], h_shelf)
  
  ## Now add variations around the "observed" thickness to form ensemble
  ## This is technically the prior

  ## Process noise parameters
  ones <- rep(1, length(domain))
  D <- rdist(domain)
  l <- 50e3
  R <- exp_cov(D, l)

  # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
  L <- t(chol(R))
  L <- as(L, "dgCMatrix")
  process_noise_info <- list(corrmat_chol = L, length_scale = l)


  df <- data.frame(x = domain, h = h)
  h_smooth <- loess(h ~ x, data = df, span = 0.05)$fitted

  h_sd <- c(rep(50, GL), rep(20, length(h) - GL))
    
  # for (i in 1:n_sims) { ## vectorise!!
  #   # h_noise <- cond_sim(mu = h_smooth$fitted, sd = h.sd, 
  #   #                                      taxis = domain, l = 50e3, nsims = 1)[, 1]
  #   h_noise <- h_sd * (process_noise_info$corrmat_chol %*% rnorm(length(domain), 0, 1))
      
  #   simulated_thickness[, i] <- h_smooth + as.vector(h_noise)
  # }

  test <- matrix(rep(h_sd, n_sims), nrow = length(h_sd), ncol = n_sims)
  Zmat <- matrix(rnorm(length(domain) * n_sims), nrow = length(domain), ncol = n_sims)
  
  h_noise <- test * (L %*% Zmat)
  h_smooth_mat <- matrix(rep(h_smooth, n_sims), nrow = length(h_smooth), ncol = n_sims)
  simulated_thickness <- h_smooth_mat + h_noise

  # for (i in 1:n_sims) {
  
  #   if (condsim_shelf) { # use conditional simulation on the ice shelf only
  #     df <- data.frame(x = domain[shelf_region], h = h[shelf_region])
  #     h.fit <- loess(h ~ x, data = df, span = 0.2)
      
  #     # h.sd <- seq(10, 1, length.out = length(domain))
  #     # shelf.sd <- h.sd[shelf_region]
  #     shelf.sd <- pmin(0.02 * h.fit$fitted, 20)
      
      
  #     simulated_shelf <- cond_sim(mu = h.fit$fitted, sd = shelf.sd, #pmin(0.02 * h.fit$fitted, 10),
  #                                 taxis = domain[shelf_region], l = 50e3, nsims = 1)
  #     simulated_thickness[, i] <- c(h[-shelf_region], simulated_shelf[, 1])
      
  #   } else { # use conditional simulation on the whole domain
  #     # Need to do this separately for grounded vs floating ice
  #     df_float <- data.frame(x = domain[shelf_region], h = h[shelf_region])
  #     h.fit_float <- loess(h ~ x, data = df_float, span = 0.2)
      
  #     df_ground <- data.frame(x = domain[-shelf_region], h = h[-shelf_region])
  #     h.fit_ground <- loess(h ~ x, data = df_ground, span = 0.05)

  #     sim_float <- cond_sim(mu = h.fit_float$fitted, sd = h.sd[shelf_region],
  #                           taxis = domain[shelf_region], l = 50e3, nsims = 1)
      
      
  #     sim_ground <- cond_sim(mu = h.fit_ground$fitted, sd = h.sd[-shelf_region],
  #                            taxis = domain[-shelf_region], l = 5e3, nsims = 1)
      
  #     simulated_thickness[, i] <- c(sim_ground[, 1], sim_float[, 1])#simulated_h[, 1]
  #     h.sd <- pmin(0.02 * c(h.fit_float$fitted, h.fit_ground$fitted), 50)
      
  #   }
  # }
  return(as.matrix(simulated_thickness))
  
}

exp_cov <- function(d, l) {
    return(exp(-3*d / l))
}