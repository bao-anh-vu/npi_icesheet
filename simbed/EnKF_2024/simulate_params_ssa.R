simulate_params <- function(n_sims, lambda, bed_peak, fixed_effect_coefs, 
                            n_bed_basis, friction_basis) {
  
  mu <- lambda$mu
  B <- lambda$B
  d <- lambda$d
  param_dim <- length(mu)
  ncol_B <- ncol(B)
  
  ## 2. Simulate parameters theta
  # theta_list <- list()
  # phi_list <- c()
  # sigma2_list <- c()
  # tau2_list <- c()
  # 
  # p1 <- proc.time()
  # z <- list()
  # eps <- list()
  # for (s in 1:n_sims) {
  #   z[[s]] <- rnorm(ncol_B, 0, 1)
  #   eps[[s]] <- rnorm(param_dim, 0, 1)
  #   theta <- mu + B %*% z[[s]] + d * eps[[s]]
  #   theta_list[[s]] <- theta
  #   phi_list[s] <- exp(theta[1]) / (1 + exp(theta[1]))
  #   sigma2_list[s] <- exp(theta[2])
  #   tau2_list[s] <- exp(theta[3])
  # 
  #   # if (n_covariates > 0) {
  #   #   beta_list[, s] <- theta[4:param_dim]
  #   # }
  # }
  # z <- unlist(z)
  # eps <- unlist(eps)
  
  p2 <- proc.time()
  z <- rnorm(ncol_B * n_sims, 0, 1)
  eps <- rnorm(param_dim * n_sims, 0, 1)
  
  # theta <- rep(mu, n_sims) + kronecker(z, B) + rep(d, n_sims) * eps
  theta <- rep(mu, n_sims) + kronecker(diag(n_sims), B) %*% z + rep(d, n_sims) * eps
  
  theta <- t(matrix(theta, nrow = param_dim, ncol = n_sims))
  theta_list <- lapply(seq_len(nrow(theta)), function(i) theta[i, ]) # split matrix rows into list
  
  ## Simulate list of beds
  
  # T matrix
  Ivec <- 0 + (domain/1000 <= domain[bed_peak]/1000) #domain[obs_locations[32]] # indicator function I(x < some location)
  Ivec2 <- 1 - Ivec
  x1 <- Ivec * (domain/1000)
  x2 <- Ivec2 * (domain/1000)
  
  alpha0 <- Ivec[bed_peak] / x2[bed_peak+1]
  alpha1 <- x1[bed_peak] / x2[bed_peak+1]
  alpha2 <- Ivec2[bed_peak+1] / x2[bed_peak+1]
  
  Tmat <- cbind(Ivec + alpha0 * x2, x1 + alpha1 * x2, Ivec2 - alpha2 * x2)
  beta <- fixed_effect_coefs
  
  # S matrix
  r <- n_bed_basis # number of basis functions
  S <- ns(domain, df = r)
  
  sim_beds <- matrix(NA, nrow = length(domain), ncol = n_sims)
  sim_friction <- matrix(NA, nrow = length(domain), ncol = n_sims)
  
  sim_sigmaK <- c()
  sim_scale <- c()
  for (s in 1:n_sims) {
    #beta <- theta_list[[s]][1:3]
    # beta12 <- (beta[1] * Ivec[bed_peak] + beta[2] * x1[bed_peak] 
    #            - beta[3] * Ivec2[bed_peak+1]) / x2[bed_peak+1]
    sigmaK <- exp(theta_list[[s]][1])
    sim_sigmaK[s] <- sigmaK
    
    eta_sim <- rnorm(r, 0, sigmaK)
    sim_beds[, s] <- Tmat %*% beta + S %*% eta_sim
    
    scale <- exp(theta_list[[s]][2])
    sim_scale[s] <- scale
    sim_friction[, s] <- scale * friction_basis
  }
  
  # sim_beds <- simulate_beds(...)
  # sim_friction <- simulate_friction(...)
  return(list(param_list = theta_list, 
              sigmaK_list = sim_sigmaK,
              scale_list = sim_scale,
              bed_list = sim_beds,
              friction_list = sim_friction))
}