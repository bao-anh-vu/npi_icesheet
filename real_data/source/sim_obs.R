sim_obs <- function(param_list, 
                    domain, 
                    phys_params, 
                    years = 20,
                    warmup = 0,
                    # steady_state, 
                    ini_velocity,
                    ini_thickness, 
                    smb = 0.5, # default of 0.5 m/yr
                    basal_melt = 0#, # default of 0 m/yr
                    # log_transform = T
                    ) { # , bed_obs) {

    ## Simulate ice thickness and velocity observations
    print("Simulating observations...")

    # if (log_transform) {
    #     param_list <- lapply(param_list, function(x) {x$friction <- exp(x$friction); x})
    # } 

    ## Measurement noise for the velocity observations ##
    ones <- rep(1, length(domain))
    D <- rdist(domain)
    l <- 5e3
    R <- exp_cov(D, l)

    # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(-sqrt(3) * D / l)
    L <- t(chol(R))
    L <- as(L, "dgCMatrix")
    msmt_noise_info <- list(corrmat_chol = L, length_scale = l)

    # sim_results <- lapply(param_list, 
    #     function(param, domain, phys_params,
    #             ini_velocity, ini_thickness, years, warmup,
    #             msmt_noise_info) {
    sim_results <- mclapply(param_list, 
        function(param, domain, phys_params,
                ini_velocity, ini_thickness, years, warmup,
                msmt_noise_info) {

        # sim_out <- solve_ssa_nl(
        #     domain = domain,
        #     bedrock = param$bedrock,
        #     friction_coef = param$friction * 1e6 * phys_params$secpera^(1 / phys_params$n),
        #     phys_params = phys_params,
        #     ini_velocity = ini_velocity,
        #     ini_thickness = ini_thickness,
        #     years = years + warmup,
        #     steps_per_yr = 100,
        #     add_process_noise = F
        # )

        sim_out <- solve_ssa_nl(
            domain = domain,
            bedrock = param$bedrock,
            friction_coef = param$friction * 1e6 * phys_params$secpera^(1 / phys_params$n),
            phys_params = phys_params,
            ini_velocity = param$ini_velocity,
            ini_thickness = param$ini_thickness,
            years = years + warmup,
            steps_per_yr = 100,
            add_process_noise = F
        )


        ## Add noise to surface to obtain surface observations
        surface_obs <- get_obs(sim_out, msmt_noise_info, warmup = warmup)

## Plot surface obs
# png("./plots/temp/surface_obs.png")
# matplot(surface_obs[,,2], type = 'l', lty = 1, col = rgb(0,0,0,0.3),)
# dev.off()
# browser()

        ## Save true thickness and velocity for comparison
        if (warmup == 0) {
            true_surface_elevs <- sim_out$all_top_surface # this is the version without added noise
            true_thicknesses <- sim_out$all_thicknesses
            true_velocities <- sim_out$all_velocities
            gl <- sim_out$grounding_line
        } else {
            true_surface_elevs <- sim_out$all_top_surface[, -(1:warmup)]
            true_thicknesses <- sim_out$all_thicknesses[, -(1:warmup)]
            true_velocities <- sim_out$all_velocities[, -(1:warmup)]
            gl <- sim_out$grounding_line[-(1:warmup)]
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
    # log_transform = log_transform,
    phys_params = phys_params,
    domain = domain,
    ini_velocity = ini_velocity,
    ini_thickness = ini_thickness,
    years = years,
    warmup = warmup,
    msmt_noise_info = msmt_noise_info,
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
        # log_transform = log_transform,
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
