initialise_ens <- function(Ne, observations, reference, data_date, 
                           # sim_thickness = TRUE, sim_bed = FALSE, sim_friction = FALSE,
                           fix_thickness = FALSE, fix_bed = FALSE, fix_friction = FALSE, 
                           condsim_shelf = TRUE, rho = 910, rho_w = 1028, 
                           # b.sill = 4000, b.nugget = 200, b.range = 50e3,
                           # alpha.sill = 8e-5, alpha.nugget = 0, alpha.range = 2.5e3, 
                           save_ini_ens = FALSE, use_basis_functions = FALSE) {
  
  domain <- reference$domain
  J <- length(domain)
  
  # 2. Simulate bed topographies
  
  print("Simulating beds...")
  
  if (use_basis_functions) {
    # simulated_beds <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/ini_bed_spline_", data_date, ".rds", sep = ""))
    
    beds <- construct_bed_basis(domain, bed_obs = observations$bed_obs,
                                obs_locations = observations$bed_obs_locations,
                                ini_params = c(-500, 0.5, 1000), n_basis = 60, 
                                parameterise_K = TRUE, use_sim_study = FALSE)
    
    simulated_beds <- matrix(rep(beds[, 1], Ne), J, Ne)
    
  } else {
    simulated_beds <- simulate_bed(nsim = Ne, 
                                   domain = domain, 
                                   obs_locations = observations$bed_obs_locations, 
                                   obs = observations$bed_obs)
  }
  
  # Plot simulated beds
  # par(mfrow = c(1,1))
  # bed_obs <- domain[observations$bed_obs_locations] / 1000
  # plot(domain/1000, simulated_beds[, 1], type = "l", xlab = "x (km)", 
  #      ylab = "Bed elevation (m)", col = "pink") # first member
  # lines(domain/1000, simulated_beds[, 2], col = "lavender") # first member
  # lines(domain/1000, simulated_beds[, 3], col = "peachpuff") # first member
  # lines(domain/1000, rowMeans(simulated_beds), lty = 1, lwd = 1.5, col = "dodgerblue") # conditional sim
  # lines(domain/1000, reference$bedrock, col = "black")
  # points(bed_obs, observations$bed_obs, col = "tomato", bg = "tomato", pch = 20)
  # legend("topright", legend = c("Mean simulated bed", "True bed"),
  #        col = c("dodgerblue", "black"), lty = 1, cex = 0.6)
  
  ## Simulate friction coefficients
  print("Simulating friction coefficient...")
  simulated_friction <- simulate_friction(nsim = Ne, domain = reference$domain) # doesn't need input, or maybe just put the mean in there
  
  # Plot simulated friction coefficients
  # par(mfrow = c(1,1))
  # plot(domain[1:500] / 1000, simulated_friction[1:500, 1], type = "l", xlab = "x (km)",
  #      ylab = "Basal friction coef (M Pa m-1/3 s1/3)", col = "pink")
  # lines(domain[1:500] / 1000, simulated_friction[1:500, 2], lty = 1, col = "lavender")
  # lines(domain[1:500] / 1000, simulated_friction[1:500, 3], lty = 1, col = "peachpuff")
  # lines(domain[1:500] / 1000, rowMeans(simulated_friction)[1:500], lty = 1, lwd = 1.5, col = "dodgerblue")
  # lines(domain[1:500] / 1000, reference$friction_coef[1:500], lty = 1, col = "black")
  # # lines(x[1:200], C[1:200], lty = 2, xlab = "x (m)",
  # # ylab = "Basal friction coef (M Pa m-1/3 s1/3)", col = "red")
  # legend("topright", legend = c("1st member", "2nd member", "Ensemble mean", "True friction"),
  #        lty = 1, col = c("lavender", "lightblue", "red", "black"), cex = 0.5)
  
  # 2.1 Transform C into alpha
  simulated_alpha <- log10(simulated_friction)
  simulated_alpha <- ifelse(is.na(simulated_alpha), 6, simulated_alpha) # replace NA values with 6
  
  # 2.2 Calculate ice thickness from z and b
  print("Calculating ice thicknesses based on simulated beds and observed top surface elevation...")
  surface_obs <- observations$surface_obs
  simulated_thickness <- matrix(0, nrow = J, ncol = Ne)
  
  for (i in 1:Ne) {
    z <- surface_obs[, 1] # use observed z at t = 0
    h <- z - simulated_beds[, i]
    
    ## Calculate GL position for each ensemble member
    GL <- gl_migrate2(z, simulated_beds[, i])
    # cat("GL: ", GL, "\n")
    
    ## Recalculate ice thickness from GL onwards (ice shelf)
    h[GL:length(h)] <- 1 / (1 - rho / rho_w) * z[GL:length(h)]
    
    if (condsim_shelf) { # use conditional simulation on the ice shelf only
      shelf_region <- GL:length(h)
      df <- data.frame(x = domain[shelf_region], h = h[shelf_region])
      h.fit <- loess(h ~ x, data = df, span = 0.2)
      
      # h.sd <- seq(10, 1, length.out = length(domain))
      # shelf.sd <- h.sd[shelf_region]
      shelf.sd <- pmin(0.02 * h.fit$fitted, 20)
      
      
      simulated_shelf <- cond_sim(mu = h.fit$fitted, sd = shelf.sd, #pmin(0.02 * h.fit$fitted, 10),
                                  taxis = domain[shelf_region], l = 50e3, nsims = 1)
      simulated_thickness[, i] <- c(h[-shelf_region], simulated_shelf[, 1])
      
    } else { # use conditional simulation on the whole domain
      
      ## Need to do this separately for grounded vs floating ice
      shelf_region <- GL:length(h)
      
      df_float <- data.frame(x = domain[shelf_region], h = h[shelf_region])
      h.fit_float <- loess(h ~ x, data = df_float, span = 0.2)

      df_ground <- data.frame(x = domain[-shelf_region], h = h[-shelf_region])
      h.fit_ground <- loess(h ~ x, data = df_ground, span = 0.05)
      
      # h.sd <- c(rep(50, GL), rep(20, length(h) - GL))
      h.sd <- pmin(0.02 * c(h.fit_float$fitted, h.fit_ground$fitted), 20)
      
      sim_float <- cond_sim(mu = h.fit_float$fitted, sd = h.sd[shelf_region],
                                  taxis = domain[shelf_region], l = 50e3, nsims = 1)
      
      
      sim_ground <- cond_sim(mu = h.fit_ground$fitted, sd = h.sd[-shelf_region],
                            taxis = domain[-shelf_region], l = 5e3, nsims = 1)
      
      simulated_thickness[, i] <- c(sim_ground[, 1], sim_float[, 1])#simulated_h[, 1]
      
    }
  }
  
  ## Assemble ice thickness, bed and friction coef into an initial ensemble
  
  # ini_ens <- matrix(0, nrow = 3*J, ncol = Ne)
  
  if (fix_thickness) { # just use reference thickness
    ini_thickness <- matrix(rep(reference$ini_thickness, Ne), J, Ne)
  } else { # use simulated thickness
    ini_thickness <- simulated_thickness
  }
  
  if (fix_bed) { #
    if (use_basis_functions) {
      ini_bed <- simulated_beds
    } else {
      ini_bed <- matrix(rep(reference$bedrock, Ne), J, Ne)  
    }
  } else {
    ini_bed <- simulated_beds
  }
  
  if (fix_friction) {
    ini_friction <- matrix(rep(log10(reference$friction_coef), Ne), J, Ne)
  } else {
    ini_friction <- simulated_alpha
  }
  
  ini_ens <- rbind(ini_thickness, ini_bed, ini_friction)
    
  # Save initial ensemble
  if (save_ini_ens) {
    saveRDS(ini_ens, file = paste("/home/babv971/SSA_model/EnKF/Output/ini_ens_", data_date, ".rds", sep = ""))
  }

  return(ini_ens)
}