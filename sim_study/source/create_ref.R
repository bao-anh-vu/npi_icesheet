##############################################################################
  #              Perturb ice sheet again to generate reference data            #
  ##############################################################################

create_ref <- function(years = 50, steps_per_yr = 26,
                        domain, bedrock, friction_coef, 
                        ini_velocity, ini_thickness, 
                        perturb_hardness = T, add_process_noise = T, seed = 123) {

  print("Perturbing ice sheet...")


#   years <- 50
#   steps_per_yr <- 26 # put these as function parameters later

#   if (!use_stored_reference) {
    t1 <- proc.time()
    ## Plug steady state velocity and ice thickness into the SSA solver again
    ## and introduce a perturbation to the ice hardness,
    ## then run for 50 years to generate "reference" data
    reference <- solve_ssa_nl(
      domain = domain,
      bedrock = bedrock,
      friction_coef = friction_coef,
      ini_velocity = ini_velocity,
      ini_thickness = ini_thickness,
      years = years, steps_per_yr = steps_per_yr,
    #   save_model_output = TRUE,
      perturb_hardness = perturb_hardness,
      add_process_noise = add_process_noise
    )
    t2 <- proc.time()

    
#   } else {
#     # Restore data form a saved file
#     # reference <- readRDS(file = "assim_data.rds")
#     reference <- readRDS(file = paste("./output/reference_", data_date, ".rds", sep = ""))
#   }


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
  par(mfrow = c(2, 2))

  # Plot bed and surface elevation -- bed not available for floating ice
  plot(x / 1000, ref_bed,
    type = "l", lwd = 2, xlim = c(0, 800), ylim = c(min(ref_bed), max(ref_top_surface)),
    xlab = "Domain (km)", ylab = "Elevation (m)"
  ) # bed elevation
  lines(x / 1000, ref_top_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # surface elevation
  lines(x / 1000, ref_bottom_surface, xlab = "Domain (km)", ylab = "Elevation (m)") # bottom surface elevation
  title("Top and bottom surface elevation")

  ## Plot GL position
  year_ind <- seq(1, years * steps_per_yr, steps_per_yr)
  plot(ref_GL_position[year_ind], 1:length(year_ind), type = "l", xlab = "Domain (km)", ylab = "Years")
  title("Grounding line position")

  ## Plot ice thickness
  thicknessRange <- c(0, max(ref_thickness, ssa_steady$current_thickness))
  fin_GL_position <- ref_GL_position[length(ref_GL_position)]
  plot_margin <- fin_GL_position - floor(fin_GL_position / 2)
  plot_region <- c(1, 800) #(fin_GL_position - plot_margin):(fin_GL_position + plot_margin) # 500:1200
  plot(x[plot_region] / 1000, ref_thickness[plot_region], type = "l", ylim = thicknessRange, xlab = "Domain (km)", ylab = "Ice thickness (m)")
  lines(x[plot_region] / 1000, ref_ini_thickness[plot_region], col = "cyan", lty = 2)
  legend("topright",
    legend = c("Current thickness", "Initial thickness"), col = c("black", "cyan"),
    lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(-0.3, -0.35)
  )
  title("Ice thickness")

  ## Plot velocity
  velocityRange <- c(0, max(max(ref_velocity), max(ssa_steady$current_velocity)))
  plot(x / 1000, ref_velocity, type = "l", ylim = velocityRange, xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  lines(x / 1000, ref_ini_velocity, col = "cyan", xlab = "Domain (km)", ylab = "Velocity (m a^-1)")
  legend("topleft",
    legend = c("Current velocity", "Initial velocity"), col = c("black", "cyan"),
    lty = 1, bty = "n", cex = 0.7, y.intersp = 0.15, inset = c(0, -0.35)
  )
  title("Horizontal velocity")

  return(reference)

}
