get_obs <- function(sim_data, msmt_noise_info, warmup = 0) {
  
    ## Velocity observations ##
    domain <- sim_data$domain
    J <- length(domain)
    gl <- sim_data$grounding_line

    if (warmup == 0) {
      all_ref_velocities <- sim_data$all_velocities # remove initial condition
      all_top_surface <- sim_data$all_top_surface # remove initial condition
    } else {
      all_ref_velocities <- sim_data$all_velocities[, -(1:warmup)] # remove initial condition
      all_top_surface <- sim_data$all_top_surface[, -(1:warmup)] # remove initial condition
    }

    obs_velocities <- all_ref_velocities # extract velocity vectors at annual resolution
    
    for (j in 1:ncol(all_ref_velocities)) {
      # vel_sd <- pmin(0.25 * all_ref_velocities[, j], 50) # Constrain stdev of measurement noise to be less than 0.25 * velocity
      # vel_sd <- 0.05 * all_ref_velocities[, j] # Constrain stdev of measurement noise to be less than 0.25 * velocity
      gl_ind <- sum(domain <= gl[j] * 1e3) # grounding line position at year j
      
      vel_sd <- c(rep(5, gl_ind), rep(50, J - gl_ind)) # set constant stdev of 5 m/s beyond grounding line
      vel_sd[vel_sd <= 0] <- 0.01 #1e-05
      
      # vel_cov <- diag(vel_sd^2, nrow(all_ref_velocities)) # measurement noise for the velocity
      
      if(sum(vel_sd < 0) > 0) {
        browser()
      }
      vel_noise <- as.numeric(vel_sd * (msmt_noise_info$corrmat_chol %*% rnorm(J, 0, 1)))
      
      # vel_noise <- rnorm(n_vel_obs, rep(0, n_vel_obs), vel_sd)
      # obs_velocities[, j] <- all_ref_velocities[, j] + mvnorm_sample(1, rep(0, nrow(all_ref_velocities)), vel_cov)
      obs_velocities[, j] <- all_ref_velocities[, j] + vel_noise #t(rmvnorm(1, rep(0, nrow(all_ref_velocities)), vel_cov))
      
    }
    
    # par(mfrow = c(2, 1))
    # plot(obs_velocities[, 2], type = "l", main = "Velocity at t = 1a")
    # plot(obs_velocities[, ncol(obs_velocities)], type = "l", main = "Velocity at t = 10a")
    # browser()

    ## 3. Top surface elevation observations ##
    # cov_surface <- diag(10^2, length(ref_top_surface)) # measurement noise for the top surface
    n_surface_obs <- nrow(all_top_surface) * ncol(all_top_surface)
    
    surface_noise <- rnorm(n_surface_obs, rep(0, n_surface_obs), rep(5, n_surface_obs))
    
    obs_surface <- all_top_surface + matrix(surface_noise, nrow = nrow(all_top_surface), ncol = ncol(all_top_surface))
    # for (j in 1:ncol(all_top_surface)) {
    #   obs_surface[, j] <- all_top_surface[, j] + t(rmvnorm(1, rep(0, nrow(all_top_surface)), cov_surface))
    #     #mvnorm_sample(1, rep(0, nrow(all_top_surface)), cov_surface)
    # }
    # obs_surface <- ref_top_surface + mvnorm_sample(1, rep(0, length(ref_top_surface)), cov_surface)
    #plot(x / 1000, surface_obs, type = "l")
    
    ## Store all observations in a list
    obs_list <- list(#  bed_obs = obs_bed,
                    #  bed_obs_locations = obs_ind,
                     surface_elev = obs_surface,#[, 2:ncol(obs_surface)],
                     velocity = obs_velocities) #[, 2:ncol(obs_velocities)])
    
    obs_arr <- abind(obs_list, along = 3)
    # if (rewrite_obs) {
    #   # saveRDS(obs_list, file = "/home/babv971/SSA_model/EnKF/Output/observations_20220307.rds")
    #   saveRDS(obs_list, file = paste("/home/babv971/SSA_model/EnKF/Output/observations_", data_date, ".rds", sep = ""))
      
    # }
    
#   } else {
#     obs_list <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/observations_", data_date, ".rds", sep = ""))
#   }
  
  # return(obs_list)
  return(obs_arr)
}