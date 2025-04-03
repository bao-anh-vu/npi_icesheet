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

create_steady_state <- function(use_stored_steady_state = FALSE,
                       rewrite_steady_state = FALSE,
                      #  use_basis_functions = FALSE,
                      #  data_date,
                       add_process_noise_in_ref = FALSE,
                       years = 1000,
                       domain,
                       friction_coef,
                       bedrock, seed = 123) {

  # if (!use_stored_steady_state) {
    t1 <- proc.time()

    C <- friction_coef # create_fric_coef(x, L) * 1e6 * (secpera)^m
    b <- bedrock

    ## Run model to steady state
    print("Running model to steady state...")

    # ssa_out <- solve_ssa_nl(domain = domain, bedrock = b, friction_coef = C, tol = 0.01, steps_per_yr = 26,
    #                         seed = ssa_seed, save_obs = TRUE, perturb_hardness = FALSE)
    ssa_out <- solve_ssa_nl(
      domain = domain, bedrock = b, friction_coef = C,
      tol = 1e-03, 
      years = years, 
      steps_per_yr = 26,
      seed = seed, 
      perturb_hardness = FALSE,
      add_process_noise = FALSE
      # add_process_noise = add_process_noise_in_ref,
      # process_noise_info = process_noise_info
    )

    ssa_steady <- ssa_out
    
    t2 <- proc.time()

   
  # } else {
  #   ## Read output from saved file
  #   ssa_steady <- readRDS(file = paste("./output/ssa_steady_", data_date, ".rds", sep = ""))
  #   # ssa_steady <- readRDS(file = "/home/babv971/SSA_model/EnKF/Output/ssa_steady_20220318.rds")
  # }

  return(ssa_steady)
  
}
