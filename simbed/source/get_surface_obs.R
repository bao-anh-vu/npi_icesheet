get_obs <- function(reference) {
  
#   if (get_new_obs) {
    
    # ## Generate synthetic observations ##
    # ref_bed <- reference$bedrock
    
    # # ref_top_surface <- reference$top_surface
    
    
    # ## 1. Bed observations ##
    
    # ## Randomly select 50 locations between x = 0 and x = 800 km
    # obs_ind <- sort(sample(length(ref_bed), n_obs))
    # obs_bed <- ref_bed[obs_ind] + rnorm(n_obs, mean = 0, sd = 20)
    
    # ## then linearly interpolate
    # plot(interpol_bed <- approx(x[obs_ind], ref_bed[obs_ind], method = "linear"),
    #      ylab = "Elevation", xlab = "x (km)", type = "l")
    # 
    # ## and add uncorrelated Gaussian noise w sigma_b = 20 m
    # cov_bed <- diag(20^2, length(interpol_bed$y))
    # bed_obs <- interpol_bed$y + mvnorm_sample(1, rep(0, length(interpol_bed$y)), cov_bed)
    # plot(interpol_bed$x / 1000, bed_obs, type = "l")
    
    ## 2. Velocity observations ##
    all_ref_velocities <- reference$all_velocities
    obs_velocities <- all_ref_velocities # extract velocity vectors at annual resolution
    n_vel_obs <- length(obs_velocities[, 1])
    
    for (j in 1:ncol(all_ref_velocities)) {
      vel_std <- pmin(0.25 * all_ref_velocities[, j], 20) # Constrain stdev of measurement noise to be less than 0.25 * velocity
      # vel_std <- pmin(0.05 * all_ref_velocities[, j], 5)
      vel_std[vel_std <= 0] <- 1e-05
      # vel_cov <- diag(vel_std^2, nrow(all_ref_velocities)) # measurement noise for the velocity
      
      if(sum(vel_std < 0) > 0) {
        browser()
      }
      
      vel_noise <- rnorm(n_vel_obs, rep(0, n_vel_obs), vel_std)
      # obs_velocities[, j] <- all_ref_velocities[, j] + mvnorm_sample(1, rep(0, nrow(all_ref_velocities)), vel_cov)
      obs_velocities[, j] <- all_ref_velocities[, j] + vel_noise #t(rmvnorm(1, rep(0, nrow(all_ref_velocities)), vel_cov))
      
      obs_velocities[1, ] <- 0 # manually set the velocity at the top to be zero
    }
    
    # par(mfrow = c(1, 2))
    # plot(obs_velocities[, 2], type = "l", main = "Velocity at t = 1a")
    # plot(obs_velocities[, ncol(obs_velocities)], type = "l", main = "Velocity at t = 50a")
    
    ## 3. Top surface elevation observations ##
    all_top_surface <- reference$all_top_surface
    # cov_surface <- diag(10^2, length(ref_top_surface)) # measurement noise for the top surface
    n_surface_obs <- nrow(all_top_surface) * ncol(all_top_surface)
    
    surface_noise <- rnorm(n_surface_obs, rep(0, n_surface_obs), rep(10, n_surface_obs))
    
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