sim_obs <- function(param_list, 
                    domain, 
                    phys_params, 
                    years = 20,
                    warmup = 0,
                    use_relaxation = F,
                    relax_years = NULL,
                    # steady_state, 
                    ini_velocity,
                    ini_surface,
                    # ini_thickness, 
                    vel_err_sd = 50, # default of 50 m/yr
                    smb = 0.5, # default of 0.5 m/yr
                    basal_melt = 0 # default of 0 m/yr
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
    #             ini_velocity, 
    #             ini_surface,
    #             ## ini_thickness, 
    #             years, warmup,
    #             use_relaxation = F,
    #             relax_years = NULL,
    #             vel_err_sd,
    #             msmt_noise_info) {

    sim_results <- mclapply(param_list, 
        function(param, domain, phys_params,
                ini_velocity, 
                ini_surface,
                ## ini_thickness, 
                years, 
                warmup,
                use_relaxation = F,
                relax_years = NULL,
                vel_err_sd,
                msmt_noise_info) {

se_grounded <- na.omit(ini_surface) # Use surface elevation in the year 2000 to initialise ice thickness
gl_ind <- length(se_grounded) # grounding line index
H_ini <- se_grounded - param$bedrock[1:gl_ind]
length_shelf <- length(domain) - length(se_grounded)

# se_gl <- se_grounded[gl_ind]
# se_shelf <- seq(from = se_gl, to = 100, length.out = length_shelf) # extend the same surface elevation at the GL to the ice shelf
# H_shelf <- - se_shelf * (params$rho_w / (params$rho_i - params$rho_w)) # thickness at grounding line based on flotation condition
# H_shelf <- se_shelf / (1 - params$rho_i / params$rho_w) # thickness at grounding line based on flotation condition
# H_shelf <- - bed_sim[(gl_ind+1):J] * params$rho_w / params$rho_i #- 100 # minus an offset to satisfy flotation condition 
H_gl <- - param$bedrock[(gl_ind+1)] * phys_params$rho_w / phys_params$rho_i #- 100 # minus an offset to satisfy flotation condition 
# H_shelf <- rep(500, length_shelf) 
H_shelf <- H_gl * exp(-0.002*seq(0, length_shelf-1))
# H_shelf <- seq(from = H_gl, to = 500, length.out = length_shelf)

# thickness_at_gl <- - bed_sim[gl_ind][1] * params$rho_w / params$rho_i
# H_shelf <- seq(thickness_at_gl - 1, 500, length.out = length_shelf)

ini_thickness <- c(H_ini, H_shelf)

# z <- get_surface_elev(H = H_ini_all, b = bed_sim, z0 = 0, rho = params$rho_i, rho_w = params$rho_w, include_GL = TRUE)



        # ini_thickness <- calculate_thickness(
        #     z = ini_surface,
        #     b = param$bedrock,
        #     z0 = 0,
        #     rho = phys_params$rho_i,
        #     rho_w = phys_params$rho_w
        # )

        sim_out <- solve_ssa_nl(
            domain = domain,
            bedrock = param$bedrock,
            friction_coef = param$friction * 1e6 * phys_params$secpera^(1 / phys_params$n),
            phys_params = phys_params,
            ini_velocity = ini_velocity,
            ini_thickness = ini_thickness,
            years = years + warmup,
            use_relaxation = use_relaxation,
            relax_years = relax_years,
            observed_thickness = ini_thickness, # this is to hold geometry "fixed" during relaxation
            steps_per_yr = 100,
            add_process_noise = F
        )

        ## Add noise to surface to obtain surface observations
        surface_obs <- get_obs(sim_out, vel_err_sd, msmt_noise_info, warmup = warmup)

# Plot surface obs
# png("./plots/temp/sim_surface_obs2.png", width = 1000, height = 600, res = 150)
# matplot(sim_out$all_top_surface[, (1:warmup)], type = 'l', lty = 1, col = "salmon",
#         # ylim = c(0, 1500),
#         xlab = "Distance along flowline (m)", ylab = "Surface elevation (m)",
#         main = "Simulated surface elevations")
# lines(sim_out$all_top_surface[, 1], col = "red")
# lines(ssa_steady$current_top_surface, col = "blue")
# matlines(sim_out$all_top_surface[, -(1:warmup)], col = "black")
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
    # ini_thickness = ini_thickness,
    ini_surface = ini_surface,
    years = years,
    warmup = warmup,
    use_relaxation = use_relaxation,
    relax_years = relax_years,
    vel_err_sd = vel_err_sd,
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
