run_sims <- function(nsims, years = 20, sim_beds = F, 
                    steady_state, bed_obs) {
    N <- nsims
    # years <- 20
    ssa_steady <- steady_state
    domain <- ssa_steady$domain

    secpera <- 31556926
    fric_scale <- 1e6 * secpera^(1 / 3)

    # 1. Simulate friction coefficients
    print("Simulating friction coefficient...")

    fric.sill <- 8e-5
    fric.nugget <- 0
    fric.range <- 10e3

    simulated_friction <- simulate_friction2(
        nsim = N, domain = domain,
        sill = fric.sill, nugget = fric.nugget,
        range = fric.range
    ) # doesn't need input, or maybe just put the mean in there

    # sim_fric_list <- lapply(1:N, function(c) simulated_friction[, c])

    ## 1.5. Simulate beds
    if (sim_beds) {
        print("Simulating beds...")
        bed_sims <- simulate_bed(
            nsim = N, domain = domain,
            obs_location = bed_obs$locations, obs = bed_obs$obs
        )

        # sim_bed_list <- lapply(1:N, function(c) bed_sims[, c])
        sim_param_list <- lapply(1:N, function(c) {
            list(
                friction = simulated_friction[, c],
                bedrock = bed_sims[, c]
            )
        })

        ## DELETE LATER: (just creating an artificial error for debugging purposes)
        # sim_param_list[[2]] <- 0
        # sim_param_list[[3]] <- 0

    } else {
        sim_param_list <- lapply(1:N, function(c) {
                    list(
                        friction = simulated_friction[, c],
                        bedrock = ssa_steady$bedrock
                    )
                })

        # sim_param_list <- lapply(1:N, function(c) list(friction = simulated_friction[, c]))
    }

    ## 2. Simulate ice thickness and velocity observations
    print("Simulating observations...")

    grounding_lines <- list()

    ## Process noise parameters
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 50e3
    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    R <- exp_cov(D, l)

    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    process_noise_info <- list(corrmat_chol = L, length_scale = l)

    t1 <- proc.time()

    # if (sim_beds) {

        # sim_results <- lapply(sim_param_list, function(sim_param) {
        sim_results <- mclapply(sim_param_list, function(sim_param) {
            
            reference <- solve_ssa_nl(
                        domain = ssa_steady$domain,
                        bedrock = ssa_steady$bedrock, 
                        friction_coef = sim_param$friction,
                        ini_velocity = ssa_steady$current_velocity,
                        ini_thickness = ssa_steady$current_thickness,
                        years = years, steps_per_yr = 52,
                        save_model_output = TRUE,
                        perturb_hardness = TRUE,
                        add_process_noise = T,
                        process_noise_info = process_noise_info
                        )
        
            ## Get surface observations (with noise added)
            surface_obs <- get_surface_obs(reference)
            
            # thickness_velocity_obs <- array(
            # surface_obs <- array(
            #     data = cbind(
            #         # reference$all_thicknesses[, 2:(years + 1)],
            #         reference$all_top_surface[, 2:(years + 1)],
            #         reference$all_velocities[, 2:(years + 1)]
            #     ),
            #     dim = c(length(domain), years, 2)
            # )

            gl <- reference$grounding_line
        # )    

            simulated_data <- list(
                # thickness_velocity_arr = thickness_velocity_obs,
                surface_obs = surface_obs,
                friction_arr = sim_param$friction,
                bed_arr = ssa_steady$bedrock,
                grounding_line = gl
            )
            # simulated_data <- obs
            return(simulated_data)
            
        },
        mc.cores = 50L,
        mc.preschedule = FALSE ## So that if one core encounters an error, the rest of the jobs run on that core will not be affected
        )

    # } else {
    #     sim_results <- mclapply(sim_param_list, function(sim_param) {
    #         # cat("sim =", sim, "\n")
    #         # steady_state = ssa_steady

    #         reference <- solve_ssa_nl(
    #             domain = ssa_steady$domain,
    #             bedrock = ssa_steady$bedrock, # use true bed
    #             friction_coef = sim_param$friction,
    #             ini_velocity = ssa_steady$current_velocity,
    #             ini_thickness = ssa_steady$current_thickness,
    #             years = years, steps_per_yr = 52,
    #             save_model_output = TRUE,
    #             perturb_hardness = TRUE,
    #             add_process_noise = T,
    #             process_noise_info = process_noise_info
    #         )

    #         thickness_velocity_obs <- array(
    #             data = cbind(
    #                 reference$all_thicknesses[, 2:(years + 1)],
    #                 reference$all_velocities[, 2:(years + 1)]
    #             ),
    #             dim = c(length(domain), years, 2)
    #         )
    #         gl <- reference$grounding_line

    #         simulated_data <- list(
    #             thickness_velocity_arr = thickness_velocity_obs,
    #             friction_arr = sim_param$friction,
    #             grounding_line = gl
    #         )
    #         # simulated_data <- obs
    #         return(simulated_data)
    #     },
    #     # steady_state = ssa_steady,
    #     mc.cores = 50L
    #     )
    # }
    
# inherits(r[[3]], "try-error")

    bad_sims <- sapply(sim_results, inherits, what = "try-error")     
    bad_sims <- which(bad_sims) # find which simulations failed
    errors <- sim_results[bad_sims]
    # errors <- str(bad_sims)

    return(list(params = sim_param_list, 
                results = sim_results, 
                bad_sims = bad_sims, 
                errors = errors))

    # t2 <- proc.time()

    # ## Note to self: should save a version with the "true" friction coef as well
    # thickness_velocity_list <- lapply(1:N, function(i) sim_results[[i]]$thickness_velocity_arr)
    # friction_list <- lapply(1:N, function(i) sim_results[[i]]$friction_arr)
    # bed_list <- lapply(1:N, function(i) sim_results[[i]]$bed_arr)
    # gl_list <- lapply(1:N, function(i) sim_results[[i]]$grounding_line)

    # concat_input <- do.call("rbind", thickness_velocity_list)
    # thickness_velocity_arr <- array(concat_input, dim = c(N, dim(thickness_velocity_list[[1]])))

    # concat_output <- do.call("rbind", friction_list)
    # friction_arr <- array(concat_output, dim = c(N, length(friction_list[[1]])))

    # if (sim_beds) {
    #     concat_output <- do.call("rbind", bed_list)
    #     bed_arr <- array(concat_output, dim = c(N, length(bed_list[[1]])))
    # }

    # concat_gl <- do.call("rbind", gl_list)
    # gl_arr <- array(concat_gl, dim = c(N, length(gl_list[[1]])))

    # if (sim_beds) {
    #     return_obj <- list(
    #         thickness_velocity_arr = thickness_velocity_arr,
    #         friction_arr = friction_arr,
    #         bed_arr = bed_arr,
    #         gl_arr = gl_arr
    #     )
    # } else {
    #     return_obj <- list(
    #         thickness_velocity_arr = thickness_velocity_arr,
    #         friction_arr = friction_arr,
    #         gl_arr = gl_arr
    #     )
    
    # }

    # return(list(input_data = return_obj,
    #             errors = errors))
}

exp_cov <- function(d, l) {
    return(exp(-3*d / l))
}
