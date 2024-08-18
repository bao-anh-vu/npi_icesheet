construct_bed_basis <- function(domain, bed_obs, obs_locations, ini_params, tol = 1e-05, 
                                n_basis = 40, n_beds = 20, parameterise_K = TRUE, 
                                use_sim_study = FALSE, use_optim = FALSE, 
                                use_poly_regr = FALSE, use_known_peak = FALSE) {
  J <- length(domain)
  nobs <- length(bed_obs)
  
  # Initial values for K_eta and beta
  beta.list <- list()
  K_eta.list <- list()
  sigma_K.list <- list()
  
  # Construct T
  
  if (use_known_peak) {
    # Ivec <- 0 + (domain/1000 <= 450)
    peak_obs <- sum(domain/1000 <= 450)
  } else if (use_poly_regr) { # use polynomial regression to find the peak
    df <- data.frame(obs_locations = domain[obs_locations], bed_elev = bed_obs)
    bed.fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.15,
                     control = loess.control(surface = "direct"))
    peak_obs <- obs_locations[which.max(bed.fit$fitted)]
    # bed_mean <- predict(bed.fit, newdata = data.frame(obs_locations = domain))
    # peak_obs <- obs_locations[as.numeric(which.max(bed_mean))] #
  } else {
    peak_obs <- obs_locations[which.max(bed_obs)]
  }
  
  Ivec <- 0 + (domain/1000 <= domain[peak_obs]/1000) #domain[obs_locations[32]] # indicator function I(x < some location)
  Ivec2 <- 1 - Ivec
  x1 <- Ivec * (domain/1000)
  x2 <- Ivec2 * (domain/1000)
  
  alpha0 <- Ivec[peak_obs] / x2[peak_obs+1]
  alpha1 <- x1[peak_obs] / x2[peak_obs+1]
  alpha2 <- Ivec2[peak_obs+1] / x2[peak_obs+1]
  
  Tmat <- cbind(Ivec + alpha0 * x2, x1 + alpha1 * x2, Ivec2 - alpha2 * x2)
  # Tmat[1:5, ]
  # dim(Tmat)
  
  ## Initial values for beta 
  beta_ini <- c(-500, 0.5, 1000) #beta_true[1:3] 
  beta12_ini <- (beta_ini[1] * Ivec[peak_obs] + beta_ini[2] * x1[peak_obs] - 
                   beta_ini[3] * Ivec2[peak_obs+1]) / x2[peak_obs+1]
  # beta_ini <- c(beta_ini, beta12_ini)
  beta.list[[1]] <- beta_ini
  
  # Construct S
  r <- n_basis # number of basis functions
  S <- ns(domain, df = r)
  par(mfrow = c(1,1))
  matplot(S[, 1:r], type = "l")
  
  sigma_K_ini <- 60
  sigma_K.list[[1]] <- sigma_K_ini 
  K_eta.list[[1]] <- diag(sigma_K_ini^2, r)
  # K_eta[[1]] <- diag(50^2, r)
  
  # Simulate some beds from these initial values
  # n_trajecs <- 50
  # b_sim_trajecs <- matrix(NA, nrow = J, ncol = n_trajecs)
  # for (k in 1:n_trajecs) {
  #   eta_sim <- t(rmvnorm(1, rep(0, r), K_eta.list[[1]]))
  #   b_sim_trajecs[, k] <- Tmat %*% beta_ini + S %*% eta_sim
  # }
  # Simulated beds
  # matplot(b_sim_trajecs, type = "l")
  # points(obs_locations, bed_obs, col = "red")
  
  # Construct C_z
  C_z <- matrix(0, nobs, J)
  ind <- cbind(1:nobs, obs_locations)
  C_z[ind] <- 1
  
  # Calculate T_z and S_z
  T_z <- C_z %*% Tmat
  S_z <- C_z %*% S
  
  # Construct sigma_epsilon
  sd_eps <- 20
  sigma_eps <- diag(sd_eps^2, nobs)
  inv_sigma_eps <- solve(sigma_eps)
  
  # plot initial bed
  # ini_bed <- Tmat %*% beta.list[[1]] + S %*% t(rmvnorm(1, rep(0, r), K_eta.list[[1]]))
  # plot(ini_bed)
  # points(obs_locations, C_z %*% ini_bed + rnorm(nobs, 0, sd_eps), col = "goldenrod")
  
  ########################### 
  ##      Run EM algo      ##
  ########################### 
  
  # iters <- 200
  l <- 2
  # curr_likelihood <- compute_log_likelihood(T_z, S_z, sigma_eps, beta[[1]], K_eta[[1]], bed_obs)
  curr_likelihood <- dmvnorm(bed_obs, mu = T_z %*% beta.list[[1]], 
                          # sigma = S_z %*% K_eta[[1]] %*% t(S_z) + sigma_eps,
                          sigma = S_z %*% K_eta.list[[1]] %*% t(S_z) + sigma_eps,
                          log = TRUE)
  # likelihood_ratio <- 10
  likelihood_diff <- 10
  
  tol <- 1e-05
  
  # for (l in 2:iters) {
  # while ( (likelihood_ratio - 1) > tol) {
  
  while(likelihood_diff > tol) {
    
    if (parameterise_K) {
      # Calculate mu_eta and sigma_eta
      # sigma_eta2 <- solve(diag(1/sigma_K.list[[l-1]]^2, r) + t(S_z) %*% inv_sigma_eps %*% S_z)
      # mu_eta2 <- sigma_eta2 %*% t(S_z) %*% inv_sigma_eps %*% (bed_obs - T_z %*% beta.list[[l-1]])
      
      K_eta <- diag(sigma_K.list[[l-1]]^2, r)
      sigma_b_inv <- solve(S_z %*% K_eta %*% t(S_z) + sigma_eps)
      sigma_eta <- K_eta - K_eta %*% t(S_z) %*% sigma_b_inv %*% S_z %*% K_eta
      mu_eta <- K_eta %*% t(S_z) %*% sigma_b_inv %*% (bed_obs - T_z %*% beta.list[[l-1]])
      
      # Update beta and K_eta
      beta.list[[l]] <- solve(t(T_z) %*% inv_sigma_eps %*% T_z) %*% t(T_z) %*% inv_sigma_eps %*% (bed_obs - S_z %*% mu_eta)
      
      if (use_optim) {
        
        opt <- optim(par = sigma_K_ini, fn = function(sigma_K, sigma_eta, mu_eta, r)
          log(det(diag(sigma_K^2, r))) +
            sum(diag(diag(1/sigma_K^2, r) %*% (sigma_eta + tcrossprod(mu_eta)))),
          method = "BFGS",
          sigma_eta = sigma_eta, mu_eta = mu_eta, r = r,
          control = list(fnscale = 1))
        
        sigma_K.list[[l]] <- opt$par
        
      } else {
        
        sigma_K.list[[l]] <- sqrt(sum(diag(sigma_eta + tcrossprod(mu_eta))) / r)  
      }
      
      K_eta.list[[l]] <- diag(sigma_K.list[[l]]^2, r)
      
    } else {
      # Calculate mu_eta and sigma_eta
      sigma_eta <- solve(solve(K_eta.list[[l-1]]) + t(S_z) %*% inv_sigma_eps %*% S_z)
      mu_eta <- sigma_eta %*% t(S_z) %*% inv_sigma_eps %*% (bed_obs - T_z %*% beta.list[[l-1]])
      
      # Update beta and K_eta
      beta_new <- solve(t(T_z) %*% inv_sigma_eps %*% T_z) %*% t(T_z) %*% inv_sigma_eps %*% (bed_obs - S_z %*% mu_eta)
      beta.list[[l]] <- beta_new
      
      K_eta.list[[l]] <- sigma_eta + mu_eta %*% t(mu_eta)
    }
    
    # Assess convergence
    # new_likelihood <- compute_log_likelihood(T_z, S_z, sigma_eps, beta[[l]], K_eta[[l]], bed_obs)
    new_likelihood <- dmvnorm(bed_obs, mu = T_z %*% beta.list[[l]], 
                           sigma = S_z %*% K_eta.list[[l]] %*% t(S_z) + sigma_eps,
                           log = TRUE)
    likelihood_diff <- new_likelihood - curr_likelihood
    # cat("Likelihood diff:", likelihood_diff, "\n")
    
    curr_likelihood <- new_likelihood
    
    l <- l+1
  }
  cat("Likelihood diff:", likelihood_diff, "\n")
  
  ## Results ## 
  beta_fin <- beta.list[[length(beta.list)]]
  beta12_fin <- (beta_fin[1] * Ivec[peak_obs] + beta_fin[2] * x1[peak_obs] 
                 - beta_fin[3] * Ivec2[peak_obs+1]) / x2[peak_obs+1]
  
  if (use_sim_study) {
    # cat("True values: beta =", beta_true, ", sigma_K =", sigma_K_true, "\n")
    df <- data.frame(ini_values = c(beta_ini, beta12_ini, sigma_K_ini), 
                     true_values = c(beta_true, sigma_K_true), 
                     estimated = c(beta_fin, 
                                   beta12_fin, sigma_K.list[[length(sigma_K.list)]]))
    rownames(df) <- c("beta01", "beta11", "beta02", "beta12", "sigma_K")
    print(df)
    
  } else {
    cat("Estimated values: beta =", c(beta_fin, beta12_fin), 
        ", sigma_K =", sigma_K.list[[length(sigma_K.list)]], "\n")
  }
  
  ## simulate beds from estimated parameters to see what they look like
  # n_beds <- 20
  bed_ens <- matrix(NA, J, ncol = n_beds)
  for (k in 1:n_beds) {
    eta_cond <- t(rmvnorm(1, mu_eta, sigma_eta)) # conditionally simulate eta using mean and variance from E-step
    bed_ens[, k] <- Tmat %*% beta.list[[length(beta.list)]] + S %*% eta_cond
  }
  
  par(mfrow = c(1,1))
  matplot(domain/1000, bed_ens, type = "l", col = "salmon", xlab = "Domain (km)", ylab = "Elevation (m)")
  lines(domain/1000, rowMeans(bed_ens), col = "red", lwd = 2)
  points(domain[obs_locations]/1000, bed_obs, col = "blue")
  
  return(list(bed_ens = bed_ens, 
              bed_peak = peak_obs,
              estimated_beta = c(beta_fin, beta12_fin),
              estimated_sigmaK = sigma_K.list[[length(sigma_K.list)]]))
  
  # saveRDS(bed_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_bed_spline_", data_date, ".rds", sep = ""))
  # saveRDS(b_true, file = paste("/home/babv971/SSA_model/EnKF/Output/true_bed_spline_", data_date, ".rds", sep = ""))
  
}