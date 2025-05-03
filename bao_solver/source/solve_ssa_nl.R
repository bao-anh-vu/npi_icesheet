## Copyright 2021 Bao Anh Vu -- not sure about the License part, might edit later
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
## http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

#' @title One-dimensional Shelfy-stream Approximation Model Solver
#' @description Numerical solution for the 1D SSA model on an ice sheet
#' @param bedrock formula identifying the dependent variable and the spatial inputs (RHS can only have one or two variables)
#' @param friction_coef basal friction coefficient (MPa m^−1/3 a^1/3) as a function of the spatial domain
#' @param velocity_bc boundary condition for the velocity (m/a) at the ice divide (x = 0)
#' @param thickness_bc boundary condition for the ice thickness (m) at the ice divide (x = 0)
#' @param tol tolerance for convergence of the ice thickness, calculated as the maximum absolute difference between the thicknesses at the previous and current iteration
#' @param years number of years for which to run the model
#' @param steps_per_yr number of timesteps per year. The timestep is calculated as dt = number of seconds per year / steps_per_yr
#' @param measure_freq frequency at which the output for surface elevation, surface velocity and ice thickness are recorded
#' @param seed (obsolete) the seed used in random midpoint displacement for generating bedrock roughness
#' @param save_model_output TRUE if saving output for surface elevation, surface velocity and ice thickness
#' @param include_GL TRUE if the grounding line should be included (the whole ice sheet will be floating if FALSE)
#' @param B_variable (obsolete) TRUE if the ice hardness is spatially variable
#' @param M_variable (obsolete) TRUE if the surface mass balance (SMB) (accumulation rate minus melt rate) is spatially variable
#' @param B_bueler (obsolete) TRUE if using the value for ice hardness in Bueler (2011)
#' @param M_bueler (obsolete) TRUE if using the value for SMB in Bueler (2011)
#' @param evolve_thickness TRUE if allowing the ice thickness to evolve in time
#' @param evolve_velocity TRUE  if allowing the velocity to evolve in time
#' @param fixed_H0 TRUE if fixing the ice thickness at x = 0 (i.e. using Dirichlet boundary condition for the ice thickness)
#' @param fixed_u0 TRUE if fixing the surface velocity at x = 0
#' @param variable_bed set to TRUE for a non-constant but piecewise linear bed topography
#' @param random_bed set to TRUE to add random noise to the bed topography
#' @param perturb_hardness TRUE if instantaneously reducing the ice hardness at t = 0 as in Gillet-Chaulet (2020)
#' @param process_noise_info list containing the process noise correlation matrix, its Cholesky factor, and the length scale
#' @return \code{solve_ssa_nl} returns an object with the following items
#' \describe{
#'  \item{"current_velocity"}{Velocity (m/a) at the final timestep }
#'  \item{"ini_velocity"}{Initial velocity profile (m/a)}
#'  \item{"all_velocities"}{Matrix containing velocity output (m/a) at specified timesteps. The output frequency can be specified using the \code{measure_freq} parameter.}
#'  \item{"current_thickness"}{Ice thickness (m) at the final timestep}
#'  \item{"ini_thickness"}{Initial ice thickness profile (m)}
#'  \item{"all_thickness"}{Matrix containing ice thickness output (m) at specified timesteps. The output frequency can be specified using the \code{measure_freq} parameter.}
#'  \item{"top_surface"}{Top surface elevation (m)}
#'  \item{"bottom_surface"}{Bottom surface elevation (m)}
#'  \item{"all_top_surface"}{Matrix containing top surface elevation output (m) at specified timesteps. The output frequency can be specified using the \code{measure_freq} parameter.}
#'  \item{"bedrock"}{Bedrock elevation (m)}
#'  \item{"friction_coef"}{Basal friction coefficient (MPa m^−1/3 a^1/3)}
#'  \item{"grounding_line"}{Vector containing grounding line positions over time} objects}
#'  \item{"domain"}{Vector of spatial mesh nodes}
#'  \item{"elapsed"}{Time elapsed while running the function}
#'  \item{"velocity_bc"}{Boundary condition for the velocity (m/a) at the ice divide (x = 0)}
#'  \item{"thickness_bc"}{Boundary condition for the ice thickness (m) at the ice divide (x = 0)}
#'  }
#' @export
#' @examples
#' C <- function (x, L) {
#'  coef <- 0.02 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) # Basal friction coefficient
#'  coef
#'  }
#' b <- rep(0, length(x))
#' \dontrun{ssa_out <- solve_ssa_nl(bedrock = b, friction_coef = C, years = 5, steps_per_yr = 26)}

################################################################################
# source("/home/babv971/SSA_model/EnKF/get_surface_elev.R")

solve_ssa_nl <- function(domain = NULL, bedrock = NULL, friction_coef = NULL,
                         velocity_bc = 0, thickness_bc = 2000,
                         ini_velocity = NULL, ini_thickness = NULL,
                         tol = 1e-06, years = 4000, steps_per_yr = 26,
                         phys_params,
                         bed_random_seed = 123,
                         save_model_output = TRUE, include_GL = TRUE,
                         evolve_thickness = TRUE, evolve_velocity = TRUE,
                         fixed_H0 = FALSE, fixed_u0 = TRUE,
                         variable_bed = TRUE, random_bed = TRUE,
                         perturb_hardness =  FALSE,
                         add_process_noise = TRUE,
                         process_noise_info = NULL) {

  t1 <- proc.time()
  
  ## Unpack physical parameters
  secpera <- phys_params$secpera
  n <- phys_params$n
  m <- phys_params$m
  rho <- phys_params$rho_i # ice density
  rho_w <- phys_params$rho_w # water density
  as <- phys_params$as # surface accumulation rate (m/s)
  ab <- phys_params$ab # melt rate (m/a)
  g <- phys_params$g
  # A <- phys_params$A # flow rate parameter
  Bg <- phys_params$B
  
  Mg <- as - ab # surface mass balance
  z0 <- 0 # ocean surface elevation

  ## Grid definition
  x <- domain
  J <- length(domain)-1 # number of grid 'spaces'
  L <- x[length(x)] # length of flowline
  dx <- x[2] - x[1] # grid spacing

  ## Basal friction coefficient

  if (!is.null(friction_coef)) {
    C <- friction_coef
  } else { # use Fabien's frictin coefficient
    C <- 0.020 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) #757.366 # Basal friction coefficient
    C <- C * 1e6 * (secpera)^m
  }

  ## Bed topography
  b <- bedrock

  ##################### Boundary conditions ########################

  u0 <- velocity_bc #ini_velocity[1] #velocity_bc # Boundary condition for the velocity at x = 0 (m/s)
  H0 <- thickness_bc #thickness_bc # BC for the ice thickness at x = 0 (m)

  ##################### Initial conditions ########################

  H_ini <- NULL
  # For the ice thickness
  if (!is.null(ini_thickness)) {
    H_ini <- ini_thickness
  } else {
    H_ini <- H0 - (H0 - 0)/(L - 0) * x
  }
  
  u_ini <- NULL
  # For the velocity
  if (!is.null(ini_velocity)) {
    u_ini <- ini_velocity #/ secpera
  } else {
    u_ini <- u0 + 0.001 / secpera * x
  }



  ###################### Finite difference ########################

  ## Predefine B(x), M(x), z(x)
  B <- c() #rep(Bg, length(x))
  M <- rep(Mg, length(x))
  z <- c()

  ##### Solve the SSA equations #####
  u_curr <- u_ini # initial guess for horizontal velocity
  H_curr <- H_ini

  ## Initial GL position
  GL <- find_gl(H_curr, b, rho_i = phys_params$rho_i, rho_w = phys_params$rho_w)
  GL_position <- GL

  # Set tolerance for fixed point iteration
  u_diff <- 999999
  H_diff <- 999999

  ## Set up output matrices
  u_mat <- matrix(0, nrow = length(u_curr), ncol = years + 1)
  H_mat <- matrix(0, nrow = length(H_curr), ncol = years + 1)
  zs_mat <- matrix(0, nrow = length(H_curr), ncol = years + 1)

  u_mat[, 1] <- as.vector(u_ini)
  H_mat[, 1] <- as.vector(H_ini)
  zs_mat[, 1] <- get_surface_elev(H_ini, b, phys_params)

  ## Iterations at which to save output
  obs_ind <- seq(0, years * steps_per_yr, steps_per_yr) # measure_freq = # times to measure per year

  tol_per_timestep <- tol / steps_per_yr # convert convergence tolerance (per annum) to per timestep
  
  i <- 1
  
  while (i <= (years * steps_per_yr) && H_diff > tol_per_timestep) {

      GL <- find_gl(H_curr, b, z0, phys_params$rho_i, phys_params$rho_w)
      
      if (i == 1 | i %% (steps_per_yr * 10) == 0) {
        cat("GL position: ", GL / J * L / 1000,  "\t")
      }

    #   if (GL <= 1) {
    #     cat("Ice sheet is no longer grounded. \n")
    #     break
    #   }
    # }

    if (evolve_velocity) {
      # Use current velocity to solve the linear system and get a new velocity
      # u_new <- solve_velocity(u_curr, H_curr, x, b, C, perturb_hardness)
      u_new <- solve_velocity(phys_params, L, J, H_curr, b, C, u_ini, velocity_bc)$u
      
      # Calculate difference between new u and current u
      u_diff <- max(abs(u_new - u_curr))

      if (i == 1 | i %% (steps_per_yr * 10) == 0) {
        cat("Change in u (m/a): ", u_diff * secpera, "\n")
      }
      # Set new velocity to current velocity
      u_curr <- u_new

      # cat("Iter", i, ": ")
      if (i == 1 | i %% (steps_per_yr * 10) == 0) {
        cat("Change in u (m/a): ", u_diff * secpera, "\t")
      }
    } else {
      # If not evolving velocity, set it to the initial velocity
      u_curr <- u_ini
    }
    
    # Now use this velocity to solve the mass balance equation
    if (evolve_thickness) {
      
      H_new <- solve_thickness(phys_params, u_curr, H_curr, x, b, 
                                steps_per_yr = steps_per_yr)
      
      if (add_process_noise && i %in% obs_ind) {
        # H_sd <- c(rep(20, GL), rep(10, length(H_new) - GL))
        H_sd <- pmin(0.02 * H_new, 20)
        H_noise <- H_sd * (process_noise_info$corrmat_chol %*% rnorm(length(x), 0, 1))
        H_new <- H_new + H_noise#$Sim1
      }

      H_diff <- max(abs(H_new - H_curr))

      png(paste0("./plots/steady_state/H_change.png"))
      plot(x/1000, H_new - H_curr, type = "l", xlab = "Domain (km)", ylab = "Ice thickness (m)")
      # lines(x/1000, H_curr, col = "salmon")
      dev.off()

      # cat("Iter", i, ": ")
      if (i == 1 | i %% (steps_per_yr * 10) == 0) {
        cat("Change in H (m/a): ", H_diff * steps_per_yr, "\t")
      }
      
      H_curr <- H_new
    }

    
    ## Store surface velocity, ice thickness and surface elevation
    if (save_model_output && i %in% obs_ind) {
      
      u_mat[, i / steps_per_yr + 1] <- as.vector(u_curr)
      H_mat[, i / steps_per_yr + 1] <- as.vector(H_curr)
      GL_position <- c(GL_position, GL) #/ J * L / 1000
      
      z <- get_surface_elev(H_curr, b, phys_params)
      zs_mat[, i / steps_per_yr + 1] <- as.vector(z)
    }
    
    ## Plot ice geometry every 100 years
    if (i == 1 | i %% (steps_per_yr * 10) == 0) {
      
      if (i != 1) {cat("Year: ", i/steps_per_yr, "\t")}

      z_curr <- get_surface_elev(H_curr, b, phys_params)
      z_b_curr <- z_curr - H_curr
      png(paste0("./plots/steady_state/z_curr", ceiling(i / steps_per_yr), ".png"))
      
      par(mfrow = c(2, 1))
      
      plot(x/1000, z_curr, type = "l", ylim = c(-1000, 2000), 
            xlab = "Domain (km)", ylab = "Elevation (m)", main = paste("Year", ceiling(i / steps_per_yr)))
      lines(x/1000, z_b_curr, col = "black")
      # lines(x/1000, zs_mat[, 1], col = "salmon") # initial ice geometry
      abline(v = x[GL]/1000, col = "black", lty = 2)
      abline(h = 0, col = "#080808", lty = 2)
      lines(x/1000, b, col = "grey")

      # png(paste0("./plots/steady_state/u_curr", ceiling(i / steps_per_yr), ".png"))
      plot(x/1000, u_curr * phys_params$secpera, type = "l", xlab = "Domain (km)", ylab = "Velocity (m/yr)")
      dev.off()

      print(head(H_curr))

    }

    # Count number of iterations
    i <- i + 1

  }
  cat("Number of years: ", floor(i/steps_per_yr), "\n")
  
  ## Calculate time taken
  t2 <- proc.time()

  ## GL position
  GL_position <- GL_position / J * L / 1000

  ## Top and bottom surface elevation
  z <- get_surface_elev(H_curr, b, z0, rho, rho_w)
  z_b <- z - H_curr

  ## Return list of output
  ssa.out <- list(current_velocity = as.vector(u_curr), ## current velocity (m/yr) # * secpera,
                  ini_velocity = u_ini, # * secpera,
                  all_velocities = u_mat, # * secpera,
                  current_thickness = H_curr,
                  ini_thickness = H_ini,
                  all_thicknesses = H_mat,
                  top_surface = z,
                  bottom_surface = z_b,
                  all_top_surface = zs_mat,
                  bedrock = b,
                  friction_coef = C,
                  grounding_line = GL_position,
                  domain = x,
                  elapsed = t2 - t1) 

  return(ssa.out)
}

