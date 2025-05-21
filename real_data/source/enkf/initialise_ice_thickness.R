initialise_ice_thickness <- function(domain, surface_obs, bed, 
                                     n_sims,
                                     rho = 910, rho_w = 1028, 
                                     condsim_shelf = FALSE, process_noise_info) {

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

  df <- data.frame(x = domain, h = h)
  h_smooth <- loess(h ~ x, data = df, span = 0.05)$fitted

  h_sd <- c(rep(50, GL), rep(20, length(h) - GL))
  h_sd_mat <- matrix(rep(h_sd, n_sims), nrow = length(h_sd), ncol = n_sims)
  Zmat <- matrix(rnorm(length(domain) * n_sims), nrow = length(domain), ncol = n_sims)
  
  L <- process_noise_info$corrmat_chol
  h_noise <- h_sd_mat * (L %*% Zmat)
  h_smooth_mat <- matrix(rep(h_smooth, n_sims), nrow = length(h_smooth), ncol = n_sims)
  simulated_thickness <- h_smooth_mat + h_noise

  return(as.matrix(simulated_thickness))
  
}

exp_cov <- function(d, l) {
    return(exp(-3*d / l))
}