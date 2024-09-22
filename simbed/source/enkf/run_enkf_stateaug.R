wendland <- function(theta, D) {
  R <- D / theta
  W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
}

run_enkf <- function(domain, years, steps_per_yr, ini_thickness, ini_bed, 
                     ini_friction_coef, ini_velocity, observations,   
                     run_analysis = TRUE, use_cov_taper = TRUE, 
                     add_process_noise = TRUE, process_noise_info,
                     transformation = "log") {
  
  print("Running filter...")
  
  Ne <- ncol(ini_thickness) # Ensemble size
  J <- length(domain)
  
  # State augmented ensemble 
  ini_ens <- rbind(ini_thickness, ini_bed, ini_friction_coef)
  enkf_means <- matrix(NA, nrow = 3*J, ncol = years + 1)
  enkf_means[, 1] <- rowMeans(ini_ens)
  enkf_ens <- list()
  enkf_ens[[1]] <- ini_ens
  ens <- ini_ens
  
  ## Velocity
  enkf_velocity_means <- matrix(0, J, years + 1)
  enkf_velocity_means[, 1] <- rowMeans(ini_velocity)
  enkf_velocities <- list()
  enkf_velocities[[1]] <- ini_velocity
  prev_velocity <- ini_velocity
  
  ## Covariance tapering
  taper_mat <- wendland(0.1 * domain[length(domain)], rdist(domain, domain))
  ens_taper1 <- kronecker(matrix(rep(1, 6), 3, 2), taper_mat)
  ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)
  
  mc.t1 <- proc.time()
  
  for (year in 2:(years+1)) {
    ##### II. Propagation #####
    cat("Forecast step", year - 1, "\n")
    # print("Propagating ensemble...")
    
    ##### II. Propagation #####
    # Append the velocities to the ensemble
    ens_all <- rbind(ens, prev_velocity)
    
    # Convert ensemble from matrix to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])
    
    # Apply the propagate function to every ensemble member
    ens.list <- lapply(ens.list, propagate, #mc.cores = 6L, 
                         domain = domain, steps_per_yr = steps_per_yr,
                         transformation = transformation)
    
    # Convert ensemble from list back to matrix
    ens_all <- matrix(unlist(ens.list), nrow = 4*J, ncol = Ne)
    
    ## Extract forecast ens and velocity ens
    # ens <- ens_all[1:(3*J), ]
    # prev_velocity <- ens_all[(3*J+1):(4*J), ]
    
    # Save the mean velocity
    # mean_prev_velocity <- rowMeans(prev_velocity)
    
    if (add_process_noise) { ## Add process noise
      # print("Adding process noise...")
      h_sd <- pmin(0.02 * rowMeans(ens_all[1:J, ]), 20)

      h_noise <- lapply(1:Ne, function(n) as.vector(h_sd * (process_noise_info$corrmat_chol %*% rnorm(length(domain), 0, 1))))
     
      h_noise_mat <- matrix(unlist(h_noise), nrow = J, ncol = Ne)  
     
      ens_all[1:J, ] <- as.matrix(ens_all[1:J, ] + h_noise_mat) # add noise to ice thickness
      
      # Convert ensemble from matrix to list
      # ens_all <- rbind(ens, prev_velocity)
      ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i]) 
    }
    
    # ens <- forecast_ens
    
    if (run_analysis) {
      ##### III. Analysis #####
      print("Running analysis...")
      
      analysis.t1 <- proc.time()
      
      # Compute P %*% t(H) and H %*% P %*% t(H)
      # HPH <- NULL #matrix(NA, nrow = 3*J, ncol = 3*J)
      # PH <- NULL #matrix(NA, nrow = 3*J, ncol = 2*J)
      
      ## Apply observation operator to every ensemble member
      HX <- mclapply(ens.list, obs_operator, transformation = transformation,
                     domain = domain, mc.cores = 6L)
      
      HX <- matrix(unlist(HX), nrow = 2*J, ncol = Ne)
      
      # Compute P %*% t(H) and H %*% P %*% t(H)
      HPH <- 1 / (Ne - 1) * tcrossprod(HX - rowMeans(HX))
      PH <- 1 / (Ne - 1) * tcrossprod(ens - rowMeans(ens), HX - rowMeans(HX))
      
      # print("Updating ensemble...")
      ## Extract yearly observations
      surf_elev_obs <- observations$surface_elev[, year]
      vel_obs <- observations$velocity[, year]
      
      y <- c(surf_elev_obs, vel_obs)

      ## Generate "observation ensemble" Y = (y, y, ..., y)
      Y <- matrix(rep(y, Ne), 2*J, Ne) 

      ## Construct measurement error covariance matrix R (depends on the velocity)
    # plot(prev_velocity) vs plot(rowMeans(HX)) here
      surf_elev_noise_sd  <- rep(10, J)
      # vel_noise_sd <- pmin(0.25 * rowMeans(prev_velocity), 20) #pmin(0.25 * vel_obs, 20)
      # vel_noise_sd[vel_noise_sd <= 0] <- 1e-05

      vel_noise_sd <- rep(20, J)
      R <- diag(c(surf_elev_noise_sd^2, vel_noise_sd^2))


      ## Calculate Kalman gain matrix K 
      # print("Calculating Kalman gain matrix...")
      if (use_cov_taper) {
        
        ##### CHANGED HERE #######
        # Tapered Kalman gain
        Sigma <- as(ens_taper2 * HPH + R, "dsCMatrix")
        # L <- t(chol(ens_taper2 * HPH + R))
        # Linv <- solve(L)
        # K <- (ens_taper1 * PH) %*% t(Linv) %*% Linv #solve(hadamard(ens_taper2, HPH) + R)
        
      } else {
        
        # K <- PH %*% solve(HPH + R) # compute the original Kalman gain matrix for comparison
        L <- t(chol(HPH + R))
        K <- PH %*% solve(t(L)) %*% solve(L)
        # Sigma <- as(HPH + R, "dsCMatrix")
        
      }
      
      ## Update ensemble
      
      # noise <- rmvn(Ne, rep(0, 2*J), R) # takes 0.4s
      noise <- matrix(rnorm(2*J*Ne, mean = 0, sd = sqrt(diag(R))), ncol = Ne) # takes 0.05s
      
      ##### CHANGED HERE #######

      if (use_cov_taper) {
        # Tapered Kalman gain
        Sigma <- as(as(as(ens_taper2 * HPH + R, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
      } else {
        # K <- PH %*% solve(HPH + R) # compute the original Kalman gain matrix for comparison
        L <- t(chol(HPH + R))
        K <- PH %*% solve(t(L)) %*% solve(L)
        # Sigma <- as(HPH + R, "dsCMatrix")
        
      }
      
      
      ens <- ens_all[1:(3*J), ] # exclude velocity from the ensemble for now
      if (use_cov_taper) {
        ens <- ens + (ens_taper1 * PH) %*% solve(a = Sigma, b = (Y - HX - noise))
      } else {
        # ens <- ens + PH %*% solve(a = Sigma, b = (Y - HX - noise))
        ens <- ens + K %*% (Y - HX - noise)
      }
      analysis.t2 <- proc.time() 
    }
    
    # Store analysis ensemble and ensemble mean
    enkf_means[, year] <- rowMeans(ens)
    enkf_ens[[year]] <- as.matrix(ens)
    
    ##### IV. Update velocity at the 100th step based on the estimated thickness #####
    ens_all <- rbind(ens, prev_velocity) # append velocity to ensemble
    
    # Convert ensemble to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])
    
    # Propagate velocity
    # might need to fix up the friction scaling in this step as it's part of the ensemble
    velocity.list <- mclapply(ens.list, 
                              function(v) { 
                                u <- as.vector(solve_velocity(prev_velocity = v[(3*J+1):(4*J)], #mean_prev_velocity,
                                                              thickness = v[1:J], 
                                                              domain = domain,
                                                              bed = v[(J+1):(2*J)],
                                                              friction = exp(v[(2*J+1):(3*J)]) * 1e6 * 31556926^(1 / 3),
                                                              perturb_hardness = TRUE))
                                
                                return(u)
                                
                              }, 
                              mc.cores = 6L)
    
    
    # Convert velocity list to matrix
    prev_velocity <- matrix(unlist(velocity.list), nrow = J, ncol = Ne)
    mean_prev_velocity <- rowMeans(prev_velocity)
    
    ## Save velocities
    enkf_velocities[[year]] <- prev_velocity
    enkf_velocity_means[, year] <- mean_prev_velocity
    
  }
  
  mc.t2 <- proc.time()
  print(mc.t2 - mc.t1)
  
  return(enkf_out = list(ens_means = enkf_means,
                         ens = enkf_ens, 
                         velocity_means = enkf_velocity_means,
                         velocities = enkf_velocities,
                         time = mc.t2 - mc.t1))  
}