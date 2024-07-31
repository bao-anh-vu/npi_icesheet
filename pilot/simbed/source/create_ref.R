# rm(list = ls())
# 
# library("Matrix")
# library("qlcMatrix")
# library("R.utils")
# 
# source("solve_ssa_nl.R")
# source("mvnorm_sample.R")
# 
# ## Seed for generating bed
# ssa_seed <- 123
# set.seed(ssa_seed)

################################################################################
# This function creates a reference ice sheet. 
# The ice sheet simulated here will be used to generate observations for
# data assimilation later.
# Input parameters for this function are currently stored inside the function,
# but this may need to change later.
# Input: currently none
# Output: reference.rds, which stores all necessary information about the 
# simulated ice sheet, including initial and current velocity, initial and curr
# ice thickness, top and bottom surface elevation, bed elevation, basal friction.
################################################################################

create_ref <- function(use_stored_steady_state = TRUE, 
                       use_stored_reference = TRUE, 
                       rewrite_steady_state = FALSE,
                       rewrite_reference = FALSE,
                       use_basis_functions = FALSE,
                       data_date,
                       add_process_noise_in_ref = FALSE,
                       friction_coef,
                       bedrock) {

  if (!use_stored_steady_state) {

t1 <- proc.time()    
    ## Parameters
    secpera <- 31556926 #seconds per annum
    n <- 3.0 # exponent in Glen's flow law
    m <- 1/n # friction exponent
    
    ## Ice shelf specifications
    L = 800e3 # length of domain
    
    ## Boundary conditions at the ice divide (x = 0)
    x0 <- 0
    u0 <- 100 / secpera # (m/s)
    H0 <- 2000 # (m)
    
    # Discretise domain
    J <- 1000 # number of steps
    dx <- L / J # increments
    x <- seq(x0, L, dx)
    
    C <- friction_coef #create_fric_coef(x, L) * 1e6 * (secpera)^m
    b <- bedrock 

    browser()
    # if (use_basis_functions) {
    #   b <- readRDS(file = paste("/home/babv971/SSA_model/EnKF/Output/true_bed_spline_", data_date, ".rds"))
    # } else {
    #   b <- create_bed(x)
    # }
    
    ## Run model to steady state
    print("Running model to steady state...")
    
    # ssa_out <- solve_ssa_nl(domain = x, bedrock = b, friction_coef = C, tol = 0.01, steps_per_yr = 26,
    #                         seed = ssa_seed, save_obs = TRUE, perturb_hardness = FALSE)
    ssa_out <- solve_ssa_nl(domain = x, bedrock = b, friction_coef = C, 
                            tol = 1e-03, steps_per_yr = 26, 
                            seed = ssa_seed, save_model_output = TRUE, 
                            perturb_hardness = FALSE, 
                            add_process_noise = T,
                            process_noise_info = process_noise_info)

    
    ssa_steady <- ssa_out
    
    if (rewrite_steady_state) {
      ## Save output
      # saveRDS(ssa_out, file = "/home/babv971/SSA_model/EnKF/Output/ssa_steady.rds")
      saveRDS(ssa_out, file = paste("./output/ssa_steady_", data_date, ".rds", sep = ""))
    }
    
  } else {
    
    ## Read output from saved file
    ssa_steady <- readRDS(file = paste("./output/ssa_steady_", data_date, ".rds", sep = ""))
    # ssa_steady <- readRDS(file = "/home/babv971/SSA_model/EnKF/Output/ssa_steady_20220318.rds")
  }

t2 <- proc.time()  
  
  browser()

  ## Extract individual components
  steady_velocity <- ssa_steady$current_velocity
  ini_velocity <- ssa_steady$ini_velocity
  steady_thickness <- ssa_steady$current_thickness
  ini_thickness <- ssa_steady$ini_thickness
  steady_top_surface <- ssa_steady$top_surface
  steady_bottom_surface <- ssa_steady$bottom_surface
  steady_bed <- ssa_steady$bedrock
  steady_friction_coef <- ssa_steady$friction_coef
  steady_GL_position <- ssa_steady$grounding_line
  x <- ssa_steady$domain
  J <- length(x)
  # all_u <- ssa_out$all_velocities
  # all_H <- ssa_out$all_thicknesses
  # all_zs <- ssa_out$all_top_surface
  
  ## Plotting
  par(mfrow = c(2,2))
  
  # Plot bed and surface elevation -- bed not available for floating ice
  plot(x/1000, steady_bed, type = "l", lwd = 2, xlim = c(0, 500), ylim = c(min(steady_bed), max(steady_top_surface)),
       xlab = "Domain (km)", ylab = "Elevation (m)") #bed elevation
  lines(x/1000, steady_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") #surface elevation
  lines(x/1000, steady_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") #bottom surface elevation
  title("Top and bottom surface elevation")
  
  ## Plot GL position
  plot(steady_GL_position, 1:length(steady_GL_position), type = "l", xlab = "Domain (km)", ylab = "Iterations")
  title("Grounding line position")
  
  ## Plot ice thickness
  fin_GL_position <- steady_GL_position[length(steady_GL_position)]
  plot_region <- (fin_GL_position - floor(J/6)):(fin_GL_position + J/2)#500:1200 ## 
  plot(x[plot_region]/1000, steady_thickness[plot_region], type = "l", ylab = "Ice thickness (m)")
  lines(x[plot_region]/1000, ini_thickness[plot_region], col = "cyan", lty = 2)
  legend("topright", legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
         lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35))
  title("Ice thickness")
  
  ## Plot velocity
  plot(x/1000, steady_velocity, type = "l", ylim = c(min(steady_velocity), max(steady_velocity)), xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  lines(x/1000, ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  legend("topleft", legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
         lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, xpd = TRUE, inset = c(0, -0.35))
  title("Horizontal velocity")
  
  ##############################################################################
  #              Perturb ice sheet again to generate reference data            #
  ##############################################################################
  print("Perturbing ice sheet...")
  
  
  years <- 50
  steps_per_yr <- 26 # put these as function parameters later
  
  if (!use_stored_reference) {
    t1 <- proc.time()
    ## Plug steady state velocity and ice thickness into the SSA solver again
    ## and introduce a perturbation to the ice hardness,
    ## then run for 50 years to generate "reference" data
    reference <- solve_ssa_nl(domain = ssa_steady$domain, 
                               bedrock = ssa_steady$bedrock, 
                               friction_coef = ssa_steady$friction_coef,
                               ini_velocity = ssa_steady$current_velocity,
                               ini_thickness = ssa_steady$current_thickness,
                               years = years, steps_per_yr = steps_per_yr,
                               seed = ssa_seed, save_model_output = TRUE, 
                               perturb_hardness = TRUE,
                               add_process_noise = add_process_noise_in_ref)
    t2 <- proc.time()
    
    if (rewrite_reference) {
      # Save reference data to a file
      saveRDS(reference, file = paste("./output/reference_", data_date, ".rds", sep = ""))
    }
    
  } else {
    # Restore data form a saved file
    # reference <- readRDS(file = "assim_data.rds")
    reference <- readRDS(file = paste("./output/reference_", data_date, ".rds", sep = ""))
    
  }
  
  
  ## Extract individual components
  ref_velocity <- reference$current_velocity
  ref_ini_velocity <- reference$ini_velocity
  ref_thickness <- reference$current_thickness
  ref_ini_thickness <- reference$ini_thickness
  ref_top_surface <- reference$top_surface
  ref_bottom_surface <- reference$bottom_surface
  ref_bed <- reference$bedrock
  ref_friction_coef <- reference$friction_coef
  ref_GL_position <- reference$grounding_line
  x <- reference$domain
  
  all_ref_velocities <- reference$all_velocities
  all_ref_thicknesses <- reference$all_thicknesses
  all_ref_top_surface <- reference$all_top_surface
  
  ## Plot reference data
  par(mfrow = c(2,2))
  
  # Plot bed and surface elevation -- bed not available for floating ice
  plot(x/1000, ref_bed, type = "l", lwd = 2, xlim = c(0, 500), ylim = c(min(ref_bed), max(ref_top_surface)),
       xlab = "Domain (km)", ylab = "Elevation (m)") #bed elevation
  lines(x/1000, ref_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") #surface elevation
  lines(x/1000, ref_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") #bottom surface elevation
  title("Top and bottom surface elevation")
  
  ## Plot GL position
  year_ind <- seq(1, years * steps_per_yr, steps_per_yr)
  plot(ref_GL_position[year_ind], 1:length(year_ind), type = "l", xlab = "Domain (km)", ylab = "Years")
  title("Grounding line position")
  
  ## Plot ice thickness
  thicknessRange <- c(0, max(ref_thickness, ssa_steady$current_thickness))
  fin_GL_position <- ref_GL_position[length(ref_GL_position)]
  plot_margin <- fin_GL_position - floor(fin_GL_position/2) 
  plot_region <- (fin_GL_position - plot_margin):(fin_GL_position + plot_margin)#500:1200
  plot(x[plot_region]/1000, ref_thickness[plot_region], type = "l", ylim = thicknessRange, xlab = "Domain (km)", ylab = "Ice thickness (m)")
  lines(x[plot_region]/1000, ref_ini_thickness[plot_region], col = "cyan", lty = 2)
  legend("topright", legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
         lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35))
  title("Ice thickness")
  
  ## Plot velocity
  velocityRange <- c(0, max(max(ref_velocity), max(ssa_steady$current_velocity)))
  plot(x/1000, ref_velocity, type = "l", ylim = velocityRange, xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  lines(x/1000, ref_ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  legend("topleft", legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
         lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(0, -0.35))
  title("Horizontal velocity")
  
  return(reference)

}

