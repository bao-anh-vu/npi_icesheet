run_enkf_missing <- function(domain, years, steps_per_yr,
                                 ini_thickness, ini_bed,
                                 ini_friction_coef, ini_velocity, observations,
                                 phys_params,
                                 correct_model_discrepancy = FALSE,
                                 model_discrepancy,
                                 missing_pattern = NULL,
                                 run_analysis = TRUE,
                                 use_cov_taper = TRUE,
                                 inflate_cov = TRUE,
                                 add_process_noise = TRUE,
                                 process_noise_info,
                                 measurement_noise_info) {
  print("Running filter...")

  # rand <- sample(c(0, 1), 1)
  #   if (rand == 1) {
  #     stop("Random error")
  #   }
  Ne <- ncol(ini_thickness) # Ensemble size
  J <- length(domain)

  mp_surface_elev <- ifelse(is.na(surf_elev_data), 0, 1)
  mp_velocity <- ifelse(is.na(velocity_data), 0, 1)

  mp_combined <- rbind(mp_surface_elev, mp_velocity)
  # C_t_ls <- lapply(1:(years+1), function(y) {
  C_t_ls <- construct_missing_matrix(missing_pattern = mp_combined)
  #   return(C_t)
  # })

  ## Store surface elevation ensemble
  se_ens_ls <- list()

  ## Model discrepancy
  se_discr_mat <- model_discrepancy$surface_elev
  vel_discr_mat <- model_discrepancy$velocity

  # State ensemble
  enkf_means <- matrix(NA, nrow = J, ncol = years + 1)
  enkf_means[, 1] <- rowMeans(ini_thickness)
  enkf_ens <- list()
  enkf_ens[[1]] <- ini_thickness
  ens <- ini_thickness

  # Parameters
  beds <- ini_bed
  friction_coefs <- ini_friction_coef

  ## Velocity
  enkf_velocity_means <- matrix(0, J, years + 1)
  enkf_velocity_means[, 1] <- rowMeans(ini_velocity)
  enkf_velocities <- list()
  enkf_velocities[[1]] <- ini_velocity
  prev_velocity <- ini_velocity

  ## Covariance tapering ##
  # taper_mat <- wendland(0.01 * domain[length(domain)], rdist(domain, domain))
  taper_mat <- wendland(0.01 * domain[length(domain)], rdist(domain, domain))

  # ens_taper1 <- kronecker(matrix(rep(1, 2), 1, 2), taper_mat)
  # ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)

  ## Likelihood ##
  llh <- c()

  mc.t1 <- proc.time()

  for (year in 1:years) {
    cat("Forecast step", year, "\n")

    ##### II. Propagation #####
    # Append the velocities to the ensemble
    ens_all <- rbind(ens, beds, friction_coefs, prev_velocity)

    # Convert ensemble from matrix to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])

    # Apply the propagate function to every ensemble member
    ens.list <- mclapply(ens.list, propagate,
      mc.cores = 10L,
      phys_params = phys_params,
      domain = domain, steps_per_yr = steps_per_yr,
      transformation = "log"
    )

    # Convert ensemble from list back to matrix
    ens_all <- matrix(unlist(ens.list), nrow = 4 * J, ncol = Ne)

    ## Extract forecast ens and velocity ens
    ens <- ens_all[1:J, ]
    prev_velocity <- ens_all[(3 * J + 1):(4 * J), ]

    #   png(paste0("plots/temp/ens_", year, ".png"))
    #   matplot(ens, type = "l")
    #   dev.off()
    # browser()

    # Save the mean velocity
    # mean_prev_velocity <- rowMeans(prev_velocity)

    if (add_process_noise) { ## Add process noise
      # print("Adding process noise...")
      h_sd <- pmin(0.02 * rowMeans(ens), 20)
      # h_sd <- pmin(0.01 * rowMeans(ens), 10)

      # h_noise <- cond_sim(mu = rep(0, J), sd = h_sd,
      #                     domain, l = 50e3, nsims = Ne)

      h_noise <- lapply(1:Ne, function(n) as.vector(h_sd * (process_noise_info$corrmat_chol %*% rnorm(length(domain), 0, 1))))

      h_noise_mat <- matrix(unlist(h_noise), nrow = J, ncol = Ne)

      ens <- as.matrix(ens + h_noise_mat) # add noise to ice thickness

      if (inflate_cov) {
        # print("Inflating covariance...")
        ens <- ens + 0.1 * (ens - rowMeans(ens))
      }

      # Convert ensemble from matrix to list
      ens_all <- rbind(ens, beds, friction_coefs, prev_velocity)
      ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])
    }

    # ## Extract yearly observations
    surf_elev_obs <- na.omit(observations$surface_elev[, year])
    vel_obs <- na.omit(observations$velocity[, year])

    y <- c(surf_elev_obs, vel_obs)

    # ens <- forecast_ens
    if (run_analysis) {
      ##### III. Analysis #####
      print("Running analysis...")

      analysis.t1 <- proc.time()

      ## Apply observation operator to every ensemble member
      HX <- mclapply(ens.list, obs_operator,
        phys_params = phys_params,
        domain = domain, transformation = "log", mc.cores = 6L
      )


      HX <- matrix(unlist(HX), nrow = 2 * J, ncol = Ne)
      HX <- C_t_ls[[year]] %*% HX # apply the missing observation matrix

      if (correct_model_discrepancy) {
        ## Add "expected" model-observation discrepancy to model output
        discr_vec <- c(se_discr_mat[, year], vel_discr_mat[, year])
        discr_vec <- C_t_ls[[year]] %*% discr_vec
        discr_ens <- matrix(rep(discr_vec, Ne), nrow = length(discr_vec), ncol = Ne)
        HX_new <- HX + discr_ens

        png(paste0("plots/temp/enkf/HX_", year, ".png"), width = 1000, height = 1000, res = 150)
        par(mfrow = c(2, 1))
        matplot(discr_ens, type = "l", col = "black")

        matplot(HX, type = "l", col = "red")
        matlines(HX_new, col = "blue")
        # plot(HX[, 1], type = "l")
        lines(y, col = "black")
        legend("bottomleft",
          legend = c("HX (red)", "HX + discrepancy", "observations"),
          col = c("red", "blue", "black"), lty = 1
        )
        dev.off()

        HX <- HX_new
      }

      # Compute P %*% t(H) and H %*% P %*% t(H)
      HPH <- 1 / (Ne - 1) * tcrossprod(HX - rowMeans(HX))
      PH <- 1 / (Ne - 1) * tcrossprod(ens - rowMeans(ens), HX - rowMeans(HX))

      # y <- y[!is.na(y)]

      ## Generate "observation ensemble" Y = (y, y, ..., y)
      Y <- matrix(rep(y, Ne), length(y), Ne)

      # ## Construct measurement error covariance matrix R (depends on the velocity)
      # surf_elev_noise_sd  <- rep(10, J)
      # vel_noise_sd <- pmin(0.25 * rowMeans(prev_velocity), 20) #pmin(0.25 * vel_obs, 20)
      # vel_noise_sd[vel_noise_sd <= 0] <- 0.1 #1e-05

      surf_elev_noise_sd <- measurement_noise_info$se_err_sd # rep(10, J)
      vel_noise_sd <- measurement_noise_info$vel_err_sd


      # vel_noise_sd <- rep(20, J)
      R <- diag(c(surf_elev_noise_sd^2, vel_noise_sd^2))
      R <- C_t_ls[[year]] %*% R %*% t(C_t_ls[[year]]) # apply the missing observation matrix

      noise <- matrix(rnorm(length(diag(R)) * Ne, mean = 0, sd = sqrt(diag(R))), ncol = Ne) # takes 0.05s

      # Calculate Kalman gain matrix K
      print("Calculating Kalman gain matrix...")
      if (use_cov_taper) {
        # Tapered Kalman gain
        ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)
        ens_taper2 <- C_t_ls[[year]] %*% ens_taper2 %*% t(C_t_ls[[year]])

        Sigma <- as(as(as(ens_taper2 * HPH + R, "dMatrix"), "symmetricMatrix"), "CsparseMatrix")
      } else {
        L <- t(chol(HPH + R))
        Linv <- solve(L)
        K <- PH %*% t(Linv) %*% Linv
      }

      ## Update ensemble
      if (use_cov_taper) {
        ens_taper1 <- kronecker(matrix(rep(1, 2), 1, 2), taper_mat)
        ens_taper1 <- ens_taper1 %*% t(C_t_ls[[year]]) ## only pick out the parts that correspond to observed data
        ens <- ens + (ens_taper1 * PH) %*% solve(a = Sigma, b = (Y - HX - noise))
      } else {
        # ens <- ens + PH %*% solve(a = Sigma, b = (Y - HX - noise))
        ens <- ens + K %*% (Y - HX - noise)
      }

      analysis.t2 <- proc.time()
    }

    png(paste0("plots/temp/ens_", year, ".png"))
    matplot(ens, type = "l")
    dev.off()


    # Store analysis ensemble and ensemble mean
    enkf_means[, year + 1] <- rowMeans(ens)
    enkf_ens[[year + 1]] <- ens

    ## Calculate surface elevation for each ensemble member
    se_ens <- apply(ens, 2, get_surface_elev, b = beds[, 1]) # all beds are the same so just use the first column

    if (correct_model_discrepancy) {
      se_discr_ens <- matrix(rep(se_discr_mat[, year], ncol(se_ens)),
        nrow = nrow(se_ens), ncol = ncol(se_ens)
      )

      se_ens_new <- se_ens + se_discr_ens

      png(paste0("plots/temp/enkf/se_ens_", year, ".png"))
      matplot(se_ens, type = "l", col = "red")
      matlines(se_ens_new, col = "blue")
      lines(observations$surface_elev[, year], col = "black")
      legend("bottomleft",
        legend = c("surface elev (red)", "surface elev + discrepancy (blue)", "observations"),
        col = c("red", "blue", "black"), lty = 1
      )
      dev.off()
      se_ens <- se_ens_new
    }

    se_ens_ls[[year + 1]] <- se_ens


    ##### IV. Update velocity at the 100th step based on the estimated thickness #####
    ens_all <- rbind(ens, beds, friction_coefs, prev_velocity) # append velocity to ensemble

    # Convert ensemble to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])

    # Propagate velocity
    velocity.list <- mclapply(
      ens.list,
      function(v) {
        # secpera <- 31556926
        # fric_scale <- 1e6 * secpera^(1 / 3)

        # evolve the velocity
        u <- as.vector(solve_velocity(
          prev_velocity = v[(3 * J + 1):(4 * J)], # mean_prev_velocity,
          thickness = v[1:J],
          phys_params = phys_params,
          domain = domain,
          bed = v[(J + 1):(2 * J)],
          friction = exp(v[(2 * J + 1):(3 * J)]) * 1e6 * 31556926^(1 / 3),
          perturb_hardness = TRUE
        ))

        return(u)
      }
    ) # ,
    # mc.cores = 6L)

    # Convert velocity list to matrix
    prev_velocity <- matrix(unlist(velocity.list), nrow = J, ncol = Ne)

    # if (sum(is.na(prev_velocity)) > 0 || sum(prev_velocity < 0) > 0) {
    #   print("Error: unrealistic velocity value")
    #   browser()
    # }
    # Correct model discrepancy
    if (correct_model_discrepancy) {
      prev_velocity <- prev_velocity + matrix(rep(vel_discr_mat[, year], ncol(prev_velocity)),
        nrow = nrow(prev_velocity), ncol = ncol(prev_velocity)
      )
    }

    mean_prev_velocity <- rowMeans(prev_velocity)

    ## Save velocities
    enkf_velocities[[year]] <- prev_velocity
    enkf_velocity_means[, year] <- mean_prev_velocity
  }

  log_likelihood <- sum(llh[-1])
  mc.t2 <- proc.time()
  print(mc.t2 - mc.t1)

  return(enkf_out = list(
    ens_means = enkf_means,
    ens = enkf_ens,
    se_ens = se_ens_ls,
    velocity_means = enkf_velocity_means,
    velocities = enkf_velocities,
    bed = beds, # [, 1], # all beds are the same anyway
    friction_coef = friction_coefs, # [, 1], # same with friction
    log_likelihood = log_likelihood,
    use_cov_taper = use_cov_taper,
    inflate_cov = inflate_cov,
    time = mc.t2 - mc.t1
  ))
}

wendland <- function(theta, D) {
  R <- D / theta
  W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
}


########################

run_enkf_missing_noH <- function(domain, years, steps_per_yr,
                                 ini_thickness, ini_bed,
                                 ini_friction_coef, ini_velocity, observations,
                                 phys_params,
                                 correct_model_discrepancy = FALSE,
                                 model_discrepancy,
                                 missing_pattern = NULL,
                                 run_analysis = TRUE,
                                 use_cov_taper = TRUE,
                                 inflate_cov = TRUE,
                                 add_process_noise = TRUE,
                                 process_noise_info,
                                 measurement_noise_info) {
  print("Running filter...")

  # rand <- sample(c(0, 1), 1)
  #   if (rand == 1) {
  #     stop("Random error")
  #   }
  Ne <- ncol(ini_thickness) # Ensemble size
  J <- length(domain)

  # if (!is.null(missing_pattern)) { # turn this part into a function later
  # mp_surface_elev <- missing_pattern$surface_elev
  # mp_velocity <- missing_pattern$vel

  ## Create missing pattern for the EnKF
  mp_surface_elev <- ifelse(is.na(surf_elev_data), 0, 1)
  mp_velocity <- ifelse(is.na(velocity_data), 0, 1)

  mp_combined <- rbind(mp_surface_elev, mp_velocity)
  # C_t_ls <- lapply(1:(years+1), function(y) {
  C_t_ls <- construct_missing_matrix(missing_pattern = mp_combined)
  #   return(C_t)
  # })

  ## Construct the "missing observation matrix" C_t here

  # mp_surface_elev_ls <- construct_missing_matrix(missing_pattern = mp_surface_elev)
  # mp_vel_ls <- construct_missing_matrix(missing_pattern = mp_velocity)

  # C_t_ls <- lapply(1:(years+1), function(y) {
  #   if (is.null(mp_surface_elev_ls[[y]])) {
  #     return(mp_vel_ls[[y]])
  #   } else if (is.null(mp_vel_ls[[y]])) {
  #     return(mp_surface_elev_ls[[y]])
  #   } else {
  #     return(as(bdiag(mp_surface_elev_ls[[y]], mp_vel_ls[[y]]), "CsparseMatrix"))
  #   }
  # )

  # png("./plots/C_t.png")
  # image(as.matrix(C_t))
  # dev.off()
  # } else {
  #   C_t_ls <- lapply(1:(years+1), function(t) as(diag(1, 2*J), "CsparseMatrix"))
  # }

  # State ensemble
  enkf_means <- matrix(NA, nrow = J, ncol = years + 1)
  enkf_means[, 1] <- rowMeans(ini_thickness)
  enkf_ens <- list()
  enkf_ens[[1]] <- ini_thickness
  ens <- ini_thickness

  # Parameters
  beds <- ini_bed
  friction_coefs <- ini_friction_coef

  ## Velocity
  enkf_velocity_means <- matrix(0, J, years + 1)
  enkf_velocity_means[, 1] <- rowMeans(ini_velocity)
  enkf_velocities <- list()
  enkf_velocities[[1]] <- ini_velocity
  prev_velocity <- ini_velocity

  ## Store surface elevation ensemble
  se_ens_ls <- list()

  ## Model discrepancy
  se_discr_mat <- model_discrepancy$surface_elev
  vel_discr_mat <- model_discrepancy$velocity

  ## Covariance tapering ##
  # taper_mat <- wendland(0.01 * domain[length(domain)], rdist(domain, domain))
  taper_mat <- wendland(0.01 * domain[length(domain)], rdist(domain, domain))

  # ens_taper1 <- kronecker(matrix(rep(1, 2), 1, 2), taper_mat)
  # ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)

  ## Likelihood ##
  llh <- c()

  mc.t1 <- proc.time()

  for (year in 1:years) {
    cat("Forecast step", year, "\n")

    ##### II. Propagation #####
    # Append the velocities to the ensemble
    ens_all <- rbind(ens, beds, friction_coefs, prev_velocity)

    # Convert ensemble from matrix to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])

    # Apply the propagate function to every ensemble member
    ens.list <- mclapply(ens.list, propagate,
      mc.cores = 10L,
      mc.preschedule = FALSE, ## So that if one core encounters an error, the rest of the jobs run on that core will not be affected
      domain = domain, steps_per_yr = steps_per_yr,
      phys_params = phys_params,
      transformation = "log"
    )

    # Convert ensemble from list back to matrix
    ens_all <- matrix(unlist(ens.list), nrow = 4 * J, ncol = Ne)

    ## Extract forecast ens and velocity ens
    ens <- ens_all[1:J, ]
    prev_velocity <- ens_all[(3 * J + 1):(4 * J), ]

    png(paste0("./plots/temp/enkf/ens_f_", year, ".png"))

    par(mfrow = c(2, 1))
    matplot(ens, type = "l")

    matplot(prev_velocity, type = "l")
    dev.off()

    # Save the mean velocity
    # mean_prev_velocity <- rowMeans(prev_velocity)

    if (add_process_noise) { ## Add process noise
      # print("Adding process noise...")
      h_sd <- pmin(0.02 * rowMeans(ens), 20)
      # h_sd <- pmin(0.01 * rowMeans(ens), 10)

      # h_noise <- cond_sim(mu = rep(0, J), sd = h_sd,
      #                     domain, l = 50e3, nsims = Ne)

      h_noise <- lapply(1:Ne, function(n) as.vector(h_sd * (process_noise_info$corrmat_chol %*% rnorm(length(domain), 0, 1))))

      h_noise_mat <- matrix(unlist(h_noise), nrow = J, ncol = Ne)

      ens <- as.matrix(ens + h_noise_mat) # add noise to ice thickness

      if (inflate_cov) {
        # print("Inflating covariance...")
        ens <- ens + 0.1 * (ens - rowMeans(ens))
      }

      # Convert ensemble from matrix to list
      ens_all <- rbind(ens, beds, friction_coefs, prev_velocity)
      ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])
    }

    ## Extract yearly observations
    surf_elev_obs <- na.omit(observations$surface_elev[, year])
    vel_obs <- na.omit(observations$velocity[, year])

    y <- c(surf_elev_obs, vel_obs)

    # ens <- forecast_ens
    if (run_analysis & length(y) > 0) { # only run analysis if observations are available
      ##### III. Analysis #####
      print("Running analysis...")

      analysis.t1 <- proc.time()

      # y <- y[!is.na(y)]

      ## Generate "observation ensemble" Y = (y, y, ..., y)
      Y <- matrix(rep(y, Ne), length(y), Ne)

      ## Construct measurement error covariance matrix R (depends on the velocity)
      # surf_elev_noise_sd  <- rep(0.5, J)
      # vel_noise_sd <- vel_err_sd #pmin(0.25 * vel_obs, 20)

      surf_elev_noise_sd <- measurement_noise_info$se_err_sd # rep(10, J)
      vel_noise_sd <- measurement_noise_info$vel_err_sd

      # vel_noise_sd <- rep(20, J)
      R <- diag(c(surf_elev_noise_sd^2, vel_noise_sd^2))
      R <- C_t_ls[[year]] %*% R %*% t(C_t_ls[[year]]) # apply the missing observation matrix

      noise <- matrix(rnorm(length(diag(R)) * Ne, mean = 0, sd = sqrt(diag(R))), ncol = Ne) # takes 0.05s

      ## Apply observation operator to every ensemble member
      HX <- mclapply(ens.list, obs_operator,
        domain = domain, phys_params = phys_params, transformation = "log", 
        mc.cores = 10L,
        mc.preschedule = FALSE ## So that if one core encounters an error, the rest of the jobs run on that core will not be affected)
      )
      HX <- matrix(unlist(HX), nrow = 2 * J, ncol = Ne)

      HX <- C_t_ls[[year]] %*% HX # apply the missing observation matrix

      if (correct_model_discrepancy) {
        ## Add "expected" model-observation discrepancy to model output
        discr_vec <- c(se_discr_mat[, year], vel_discr_mat[, year])
        discr_vec <- C_t_ls[[year]] %*% discr_vec
        discr_ens <- matrix(rep(discr_vec, Ne), nrow = length(discr_vec), ncol = Ne)
        HX_new <- HX + discr_ens

        png(paste0("plots/temp/enkf/HX_", year, ".png"), width = 1000, height = 1000, res = 150)
        par(mfrow = c(2, 1))
        matplot(discr_ens, type = "l", col = "black")

        matplot(HX, type = "l", col = "red")
        matlines(HX_new, col = "blue")
        # plot(HX[, 1], type = "l")
        lines(y, col = "black")
        legend("bottomleft",
          legend = c("HX (red)", "HX + discrepancy", "observations"),
          col = c("red", "blue", "black"), lty = 1
        )
        dev.off()

        HX <- HX_new
      }

      HA <- HX - 1 / Ne * (HX %*% rep(1, Ne)) %*% t(rep(1, Ne))

      Pa <- 1 / (Ne - 1) * HA %*% t(HA) + R ## this bit causes mclapply to struggle

      L_Pa <- t(chol(Pa))
      Linv <- solve(L_Pa)

      if (use_cov_taper) {
        ens_taper2 <- kronecker(matrix(rep(1, 4), 2, 2), taper_mat)
        ens_taper2 <- C_t_ls[[year]] %*% ens_taper2 %*% t(C_t_ls[[year]])
        Pa <- ens_taper2 * Pa
      }

      A <- ens - 1 / Ne * (ens %*% rep(1, Ne)) %*% t(rep(1, Ne))
      analysis_update <- 1 / (Ne - 1) * A %*% t(HA) %*% t(Linv) %*% Linv %*% (Y - HX - noise)
      ens_a <- ens + analysis_update

      png(paste0("plots/temp/enkf/ens_a_", year, ".png"))
      matplot(ens, type = "l")
      matlines(ens_a, col = "red")
      dev.off()
      ens <- ens_a
      analysis.t2 <- proc.time()
    }



    ## Calculate surface elevation for each ensemble member
    se_ens <- apply(ens, 2, get_surface_elev, b = beds[, 1]) # all beds are the same so just use the first column

    if (correct_model_discrepancy) {
      se_discr_ens <- matrix(rep(se_discr_mat[, year], ncol(se_ens)),
        nrow = nrow(se_ens), ncol = ncol(se_ens)
      )

      se_ens_new <- se_ens + se_discr_ens

      png(paste0("plots/temp/enkf/se_ens_", year, ".png"))
      matplot(se_ens, type = "l", col = "red")
      matlines(se_ens_new, col = "blue")
      lines(observations$surface_elev[, year], col = "black")
      legend("bottomleft",
        legend = c("surface elev (red)", "surface elev + discrepancy (blue)", "observations"),
        col = c("red", "blue", "black"), lty = 1
      )
      dev.off()
      se_ens <- se_ens_new
    }

    se_ens_ls[[year + 1]] <- se_ens

    # Store analysis ensemble and ensemble mean
    enkf_means[, year + 1] <- rowMeans(ens)
    enkf_ens[[year + 1]] <- ens

    ##### IV. Update velocity at the 100th step based on the estimated thickness #####
    ens_all <- rbind(ens, beds, friction_coefs, prev_velocity) # append velocity to ensemble

    # Convert ensemble to list
    ens.list <- lapply(seq_len(ncol(ens_all)), function(i) ens_all[, i])



    # Propagate velocity
    velocity.list <- mclapply(ens.list,
      function(v) {
        # secpera <- 31556926
        # fric_scale <- 1e6 * secpera^(1 / 3)

        # evolve the velocity
        u <- as.vector(solve_velocity(
          prev_velocity = v[(3 * J + 1):(4 * J)], # mean_prev_velocity,
          thickness = v[1:J],
          domain = domain,
          bed = v[(J + 1):(2 * J)],
          friction = exp(v[(2 * J + 1):(3 * J)]) * 1e6 * 31556926^(1 / 3),
          phys_params = phys_params
        ))

        return(u)
      },
      mc.cores = 10L,
      mc.preschedule = FALSE ## So that if one core encounters an error, the rest of the jobs run on that core will not be affected
    )

    # Convert velocity list to matrix
    prev_velocity <- matrix(unlist(velocity.list), nrow = J, ncol = Ne)

    # Correct model discrepancy
    if (correct_model_discrepancy) {
      prev_velocity <- prev_velocity + matrix(rep(vel_discr_mat[, year], ncol(prev_velocity)),
        nrow = nrow(prev_velocity), ncol = ncol(prev_velocity)
      )
    }

    mean_prev_velocity <- rowMeans(prev_velocity)

    png(paste0("plots/temp/enkf/vel_ens_", year, ".png"))
    matplot(prev_velocity, type = "l")
    lines(observations$velocities[, year], col = "red")
    dev.off()

    # if (sum(is.na(prev_velocity)) > 0 || sum(prev_velocity < 0) > 0) {
    #   print("Error: unrealistic velocity value")
    #   browser()
    # }

    ## Save velocities
    enkf_velocities[[year + 1]] <- prev_velocity
    enkf_velocity_means[, year + 1] <- mean_prev_velocity
  }

  log_likelihood <- sum(llh[-1])
  mc.t2 <- proc.time()
  print(mc.t2 - mc.t1)

  return(enkf_out = list(
    ens_means = enkf_means,
    ens = enkf_ens,
    surf_elev_ens = se_ens_ls,
    velocity_means = enkf_velocity_means,
    velocities = enkf_velocities,
    bed = beds, # [, 1], # all beds are the same anyway
    friction_coef = friction_coefs, # [, 1], # same with friction
    log_likelihood = log_likelihood,
    use_cov_taper = use_cov_taper,
    inflate_cov = inflate_cov,
    time = mc.t2 - mc.t1
  ))
}

# wendland <- function(theta, D) {
#   R <- D / theta
#   W <- (R <= 1) * (1 - R)^4 * (1 + 4 * R)
# }
