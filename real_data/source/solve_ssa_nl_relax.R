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
#' @param seed the seed used in random midpoint displacement for generating bedrock roughness
#' @param include_GL TRUE if the grounding line should be included (the whole ice sheet will be floating if FALSE)
#' @param evolve_thickness TRUE if allowing the ice thickness to evolve in time
#' @param evolve_velocity TRUE  if allowing the velocity to evolve in time
#' @param fixed_H0 TRUE if fixing the ice thickness at x = 0 (i.e. using Dirichlet boundary condition for the ice thickness)
#' @param fixed_u0 TRUE if fixing the surface velocity at x = 0
#' @param variable_bed set to TRUE for a non-constant but piecewise linear bed topography
#' @param random_bed set to TRUE to add random noise to the bed topography
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
#'  coef <- 0.025 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) # Basal friction coefficient
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
                        #  measure_freq = 1, bed_random_seed = 123,
                         include_GL = TRUE,
                         evolve_thickness = TRUE, evolve_velocity = TRUE,
                         fixed_H0 = FALSE, fixed_u0 = TRUE,
                        #  variable_bed = TRUE, random_bed = TRUE,
                         add_process_noise = FALSE,
                         process_noise_info = NULL,
                         use_relaxation = FALSE,
                         observed_thickness = NULL,
                         relax_rate = 0
                        #  basal_melt = 0, smb = 0.5
                         ) {

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

  # ## Define physical parameters
  # secpera <- 31556926 #seconds per annum
  # n <- 3.0 # exponent in Glen's flow law
  # m <- 1/n # friction exponent
  # rho <- 910.0 #ice density
  # rho_w <- 1028.0 #sea water density
  # #r <- rho / rho_w # define this ratio for convenience
  # g <- 9.81 #gravity constant
  # A <- 1.4579e-25 #ice hardness
  # z0 <- 0 # ocean surface elevation


  ############################### Domain ###############################
  x <- domain
  L <- x[length(x)]
  # dx <- x[2] - x[1]
  dx <- mean(x[2:length(x)] - x[1:(length(x)-1)]) # I had to do this because there were some small discrepancies in the spacing when I regridded the flowline
  J <- round(L / dx)

  ## Basal friction coefficient

  if (!is.null(friction_coef)) {
    C <- friction_coef
  } else {
    C <- 0.020 + 0.015 * sin(5 * 2 * pi * x / L) * sin(100 * 2 * pi * x / L) #757.366 # Basal friction coefficient
    C <- C * 1e6 * (secpera)^m
  }

  ## Bed topography
  b <- rep(0, length(x))
  if (!is.null(bedrock)) {
    b <- bedrock
  } else {

    if (variable_bed) {
      left <- (x <= 450e3)
      right <- (x > 450e3)
      b[left] <- - 1100 + x[left] / 1000
      b[right] <- - 650 - 5 * (x[right] / 1000 - 450)
      b <- b + 500
    } else {
      b <- rep(-650, length(x))
    }

    if (random_bed) { # Use random midpoint displacement to generate bed roughness
      set.seed(bed_random_seed)

      # Parameters for random midpoint displacement
      K <- 2
      br <- rep(0, K)
      mu <- 0
      sigma <- 500
      h <- 0.7
      reps <- 12

      for (i in 1:reps) {
        midpoints <- 0.5 * (br[2:length(br)] + br[1:(length(br) - 1)])
        midpoints <- midpoints + rnorm(length(midpoints), mu, sigma)
        indices <- seq(2, length(br), 1)
        br <- insert(br, indices, midpoints)
        sigma <- sigma / (2^h)
      }
      b <- b + br[1:length(b)]
    }
  }

  ##################### Boundary conditions ########################

  u0 <- velocity_bc #ini_velocity[1] #velocity_bc # Boundary condition for the velocity at x = 0 (m/s)
  H0 <- thickness_bc #thickness_bc # BC for the ice thickness at x = 0 (m)

  ##################### Initial conditions ########################

  H_ini <- NULL
  # For the ice thickness
  if (!is.null(ini_thickness)) {
    H_ini <- ini_thickness
  } else {
    elev_at_tail <- 500 
    H_ini <- H0 - (H0 - elev_at_tail)/(L - 0) * x
    # H_ini <- H0 - (H0 - elev_at_tail)/L^2 * x^2
  }

  u_ini <- NULL
  # For the velocity
  if (!is.null(ini_velocity)) {
    u_ini <- ini_velocity #/ secpera
  } else {
    u_ini <- u0 + 0.001 * x
  }

  ###################### Finite difference ########################

  ## Predefine B(x), M(x), z(x)
  B <- c() #rep(Bg, length(x))
  # M <- rep(Mg, length(x))
  z <- c()

  ## The equation we're trying to solve is of the form
  ## W(x) du/dx + alpha(x) u = beta(x)
  ## Set up some vectors for alpha, beta, W

  alpha <- rep(0, length(x))
  W <- rep(0, length(x))
  beta <- rep(0, length(x))

  ## Define a strain rate for regularisation of du/dx
  eps_reg <- (1.0 / secpera) / L # strain rate of 1 m/a over length of shelf


  ## Implement fixed point iteration

  u_curr <- u_ini # initial guess for horizontal velocity
  H_curr <- H_ini
  GL_position <- c()

  # Set tolerance for fixed point iteration
  u_diff <- 999999
  H_diff <- 999999

  ## Set up output matrices
  u_mat <- matrix(0, nrow = length(u_curr), ncol = years + 1)
  H_mat <- matrix(0, nrow = length(H_curr), ncol = years + 1)
  zs_mat <- matrix(0, nrow = length(H_curr), ncol = years + 1)

  u_mat[, 1] <- as.vector(u_ini)
  H_mat[, 1] <- as.vector(H_ini)
  zs_mat[, 1] <- get_surface_elev(H = H_ini, b = b, rho = rho, rho_w = rho_w)

  ## Years to save output
  obs_ind <- seq(0, years * steps_per_yr, steps_per_yr) # measure_freq = # times to measure per year

  tol_per_timestep <- tol / steps_per_yr # convert convergence tolerance (per annum) to per timestep
  
  i <- 1
  
  while (i <= (years * steps_per_yr) && H_diff > tol_per_timestep) {
    
    if (include_GL) {
      GL <- gl_migrate(H_curr, b, z0, rho, rho_w)
      # GL_position <- c(GL_position, GL) #/ J * L / 1000
      
      if (GL <= 1) {
        cat("Ice sheet is no longer grounded. \n")
        break
      }
    }
    
    ## Plot ice geometry every 100 years
    if (i == 1 | i %% (steps_per_yr * 10) == 0) {
    # if (i == 1 | i %% 100 == 0) {
      
      cat("Year: ", i/steps_per_yr, "\t")
      
      z_curr <- get_surface_elev(H_curr, b, z0, rho, rho_w)
      z_b_curr <- z_curr - H_curr

cat("GL position: ", GL / J * L / 1000,  "\t")
      
if (use_relaxation) {
  plot_name <- "z_curr_relax"
} else {
  plot_name <- "z_curr"
}

      png(paste0("./plots/steady_state/", plot_name, ceiling(i / steps_per_yr), ".png"))
      # png(paste0("./plots/steady_state/z_curr", i, ".png"))
      
      par(mfrow = c(2, 1))
      
      plot(x/1000, z_curr, type = "l", ylim = c(-2000, 2000), 
            xlab = "Domain (km)", ylab = "Elevation (m)", main = paste("Year", ceiling(i / steps_per_yr)))
      lines(x/1000, z_b_curr, col = "black")
      lines(x/1000, zs_mat[, 1], col = "salmon") # initial ice geometry
      abline(v = x[GL]/1000, col = "black", lty = 2)
      abline(v = x[GL_position[1]]/1000, col = "salmon", lty = 2) # initial GL
      abline(h = 0, col = "turquoise", lty = 2)
      lines(x/1000, b, col = "grey")

      # png(paste0("./plots/steady_state/u_curr", ceiling(i / steps_per_yr), ".png"))
      plot(x/1000, u_curr, type = "l", xlab = "Domain (km)", ylab = "Velocity (m/yr)")
      lines(x/1000, u_mat[, 1], col = "salmon") # initial velocity
      abline(v = x[GL]/1000, col = "black", lty = 2)
      abline(v = x[GL_position[1]]/1000, col = "salmon", lty = 2) # initial GL
      legend("bottomright", legend = c("Current", "Initial"), col = c("black", "salmon"), lty = 1)
      dev.off()

    }

    # Now use this velocity to solve the mass balance equation
    if (evolve_thickness) {

      if (use_relaxation) {
        # relax_years <- 200 #00
        tau <- 50 # relaxation timescale in years
        relax_rate <- (observed_thickness - H_curr) / tau # m/yr
      }

      H_new <- solve_thickness(u_curr, H_curr, x, b, steps_per_yr = steps_per_yr,
                                # use_relaxation = use_relaxation, 
                                # relax_thickness = H_ini, 
                                relax_rate = relax_rate,
                                as = as, ab = ab)

      if (add_process_noise && i %in% obs_ind) {
        H_sd <- pmin(0.02 * H_new, 20)
        
        H_noise <- H_sd * (process_noise_info$corrmat_chol %*% rnorm(length(x), 0, 1))
        
        H_new <- H_new + H_noise#$Sim1
        
      }
      
      H_diff <- max(abs(H_new - H_curr))

      # browser()
      # cat("Iter", i, ": ")
      # if (i %% 1000 == 0) {
        # cat("Change in H (m): ", H_diff, "\n")
      # }
      
      H_curr <- H_new
    }

    if (evolve_velocity) {
      # Use current velocity to solve the linear system and get a new velocity
      # u_new <- solve_velocity_og(u_curr, H_curr)
      u_new <- solve_velocity(u_curr, H_curr, x, b, C, phys_params = params)
      
      # Calculate difference between new u and old u
      u_diff <- max(abs(u_new - u_curr))
      # cat("Change in u (m/a): ", u_diff * secpera, "\n")

      # Set new velocity to current velocity
      u_curr <- u_new

# browser()

    }
    
    
    ## Store surface velocity, ice thickness and surface elevation
    if (i %in% obs_ind) {
      # u_mat[, ceiling(i / steps_per_yr)] <- as.vector(u_curr)
      # H_mat[, ceiling(i / steps_per_yr)] <- as.vector(H_curr)
      
      u_mat[, i / steps_per_yr + 1] <- as.vector(u_curr)
      H_mat[, i / steps_per_yr + 1] <- as.vector(H_curr)
      GL_position <- c(GL_position, GL) #/ J * L / 1000

      z <- get_surface_elev(H_curr, b, z0, rho, rho_w)
      # zs_mat[, ceiling(i / steps_per_yr)] <- as.vector(z)
      zs_mat[, i / steps_per_yr + 1] <- as.vector(z)
    }
    
    # Count number of iterations
    i <- i + 1

    # if (!is.null(iters)) {
    #   stop_condition <- i < iters
    # } else {
    #   stop_condition <- u_diff > tol
    # }
  }
  # cat("Number of iterations: ", i - 1, "\n")
  
  ## Calculate time taken
  t2 <- proc.time()

  ## GL position
  GL_position <- GL_position / J * L / 1000

  ## Top and bottom surface elevation
  z <- get_surface_elev(H_curr, b, z0, rho, rho_w)
  z_b <- z - H_curr

  ## Truncate matrices of velocity, ice thickness and surface elevation to the number of years it takes for convergence
  u_mat <- u_mat[, 1:(floor(i / steps_per_yr + 1) )]
  H_mat <- H_mat[, 1:(floor(i / steps_per_yr + 1) )]
  zs_mat <- zs_mat[, 1:(floor(i / steps_per_yr + 1))]

  ## Return list of output
  ssa.out <- list(current_velocity = as.vector(u_curr), ## current velocity (m/yr) # * secpera,
                  ini_velocity = u_ini, # * secpera,
                  all_velocities = u_mat, # * secpera,
                  current_thickness = H_curr,
                  ini_thickness = H_ini,
                  all_thicknesses = H_mat,
                  current_top_surface = z,
                  current_bottom_surface = z_b,
                  all_top_surface = zs_mat,
                  bedrock = b,
                  friction_coef = C,
                  grounding_line = GL_position,
                  domain = x,
                  elapsed = t2 - t1) ## might want to output some of the input as well, e.g. initial height

  ssa.out
}




