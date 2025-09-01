sim_obs <- function(param_list, years = 20,
                    warmup = 0,
                    # steady_state, 
                    ini_velocity,
                    ini_thickness, 
                    smb = 0.5, # default of 0.5 m/yr
                    basal_melt = 0, # default of 0 m/yr
                    log_transform = T) { # , bed_obs) {
    # N <- nsims
    # # years <- 20
    # ssa_steady <- steady_state
    # domain <- ssa_steady$domain

    # # 1. Simulate friction coefficients
    # print("Simulating friction coefficient...")

    # fric.sill <- 8e-5
    # fric.nugget <- 0
    # fric.range <- 10e3

    # simulated_friction <- simulate_friction2(
    #     nsim = N, domain = domain,
    #     sill = fric.sill, nugget = fric.nugget,
    #     range = fric.range
    # ) # doesn't need input, or maybe just put the mean in there

    # # sim_fric_list <- lapply(1:N, function(c) simulated_friction[, c])

    # ## 1.5. Simulate beds
    # if (sim_beds) {
    #     print("Simulating beds...")

    #     bed_sims <- simulate_bed(
    #         nsim = N, domain = domain,
    #         obs_location = bed_obs$locations, obs = bed_obs$obs
    #     )

    #     # sim_bed_list <- lapply(1:N, function(c) bed_sims[, c])
    #     param_list <- lapply(1:N, function(c) {
    #         list(
    #             friction = simulated_friction[, c],
    #             bedrock = bed_sims[, c]
    #         )
    #     })

    #     ## DELETE LATER: (just creating an artificial error for debugging purposes)
    #     # param_list[[2]] <- 0
    #     # param_list[[3]] <- 0

    # } else {
    #     param_list <- lapply(1:N, function(c) list(friction = simulated_friction[, c]))
    # }


    ## 2. Simulate ice thickness and velocity observations
    print("Simulating observations...")

    # grounding_lines <- list()

    ## Process noise parameters (for ice thickness)
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 50e3
    R <- exp_cov(D, l)

    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    process_noise_info <- list(corrmat_chol = L, length_scale = l)

    # print("Simulating ground truth...")
    # reference <- create_ref(use_stored_steady_state = F,
    #                         use_stored_reference = F,
    #                         add_process_noise_in_ref = T,
    #                         rewrite_steady_state = F,
    #                         rewrite_reference = F,
    #                         data_date = data_date,
    #                         friction_coef = simulated_friction[, 1],
    #                         bedrock = bed_sims[, 1]
    #                         )

    ## Change friction unit to M Pa m^1/3 s^-1/3
    secpera <- 31556926
    fric_scale <- 1e6 * secpera^(1 / 3)

    # sim_results <- lapply(param_list, function(param, log_transform) {
    sim_results <- mclapply(param_list, function(param, log_transform) {

        # print("Simulating data...")
        ## Artificial error for debugging purposes
        # s <- sample(1:2)
        # stopifnot(s == 1)

        if (log_transform) {
            fric <- exp(param$friction)
        } else {
            fric <- param$friction
        }

        reference <- solve_ssa_nl(
            domain = ssa_steady$domain,
            bedrock = param$bedrock,
            friction_coef = fric * fric_scale,
            # ini_velocity = ssa_steady$current_velocity,
            # ini_thickness = ssa_steady$current_thickness,
            ini_velocity = ini_velocity,
            ini_thickness = ini_thickness,
            years = years, steps_per_yr = 100,
            # save_model_output = TRUE,
            # perturb_hardness = FALSE,
            add_process_noise = T,
            process_noise_info = process_noise_info,
            smb = smb_avg,
            basal_melt = basal_melt
        )

        ## Get surface observations (with noise added)
        surface_obs <- get_obs(reference)

        ## Save true thickness and velocity for comparison
        if (warmup == 0) {
            true_surface_elevs <- reference$all_top_surface
            true_thicknesses <- reference$all_thicknesses
            true_velocities <- reference$all_velocities
            gl <- reference$grounding_line
        } else {
            true_surface_elevs <- reference$all_top_surface[, -(1:warmup)]
            true_thicknesses <- reference$all_thicknesses[, -(1:warmup)]
            true_velocities <- reference$all_velocities[, -(1:warmup)]
            gl <- reference$grounding_line[-(1:warmup)]
        }

        simulated_data <- list(
            true_surface_elevs = true_surface_elevs,
            true_thicknesses = true_thicknesses,
            true_velocities = true_velocities,
            surface_obs = surface_obs,
            friction_arr = param$friction,
            bed_arr = param$bedrock,
            grounding_line = gl
        )
        # simulated_data <- obs
        return(simulated_data)
    },
    log_transform = log_transform,
    mc.cores = 50L,
    ## mc.allow.fatal = TRUE,
    mc.preschedule = FALSE ## So that if one core encounters an error, the rest of the jobs run on that core will not be affected
    )

    # inherits(r[[3]], "try-error")

    bad_sims <- sapply(sim_results, inherits, what = "try-error")
    bad_sims <- which(bad_sims) # find which simulations failed
    errors <- sim_results[bad_sims]
    # errors <- str(bad_sims)

    return(list(
        params = param_list,
        log_transform = log_transform,
        results = sim_results,
        bad_sims = bad_sims,
        errors = errors
    ))

    # t2 <- proc.time()


    # return(list(input_data = return_obj,
    #             errors = errors))
}

exp_cov <- function(d, l) {
    return(exp(-3 * d / l))
}
