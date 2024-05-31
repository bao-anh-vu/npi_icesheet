run_sims <- function(nsims) {

    N <- nsims
    years <- 50
    domain <- ssa_steady$domain

    secpera <- 31556926
    fric_scale <- 1e6 * secpera^(1/3)

    # 1. Simulate friction coefficients
    print("Simulating friction coefficient...")

    fric.sill <- 8e-5
    fric.nugget <- 0
    fric.range <- 10e3

    simulated_friction <- simulate_friction2(nsim = N, domain = domain,
                                            sill = fric.sill, nugget = fric.nugget, 
                                            range = fric.range) # doesn't need input, or maybe just put the mean in there

    sim_fric_list <- lapply(1:N, function(c) simulated_friction[, c])

    ## 2. Simulate ice thickness and velocity observations
    print("Simulating observations...")

    grounding_lines <- list()

    ## Process noise parameters
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 50e3
    R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    process_noise_info <- list(corrmat_chol = L, length_scale = l)

    t1 <- proc.time()
    sim_results <- mclapply(sim_fric_list, function(simulated_friction) {
    # cat("sim =", sim, "\n")
    # steady_state = ssa_steady
    reference <- solve_ssa_nl(domain = ssa_steady$domain, 
                                bedrock = ssa_steady$bedrock, 
                                friction_coef = simulated_friction,
                                ini_velocity = ssa_steady$current_velocity,
                                ini_thickness = ssa_steady$current_thickness,
                                years = years, steps_per_yr = 52, 
                                save_model_output = TRUE, 
                                perturb_hardness = TRUE,
                                add_process_noise = T,
                                process_noise_info = process_noise_info)    

    thickness_velocity_obs <- array(data = cbind(reference$all_thicknesses[, 2:(years+1)], 
                                reference$all_velocities[, 2:(years+1)]),
                    dim = c(length(domain), years, 2))
    gl <- reference$grounding_line

    simulated_data <- list(thickness_velocity_arr = thickness_velocity_obs,
                            friction_arr = simulated_friction,
                            grounding_line = gl)
    # simulated_data <- obs
    return(simulated_data)
    },
    #steady_state = ssa_steady, 
    mc.cores = 50L)
    t2 <- proc.time()

    ## Note to self: should save a version with the "true" friction coef as well
    thickness_velocity_list <- lapply(1:N, function(i) sim_results[[i]]$thickness_velocity_arr)
    friction_list <- lapply(1:N, function(i) sim_results[[i]]$friction_arr)
    gl_list <- lapply(1:N, function(i) sim_results[[i]]$grounding_line)

    concat_input <- do.call("rbind", thickness_velocity_list)
    thickness_velocity_arr <- array(concat_input, dim = c(N, dim(thickness_velocity_list[[1]])))

    concat_output <- do.call("rbind", friction_list)
    friction_arr <- array(concat_output, dim = c(N, length(friction_list[[1]]), 1L, 1L))

    concat_gl <- do.call("rbind", gl_list)
    gl_arr <- array(concat_gl, dim = c(N, 1L, length(gl_list[[1]]), 1L))

    return(list(thickness_velocity_arr = thickness_velocity_arr, 
                friction_arr = friction_arr, 
                gl_arr = gl_arr))
}