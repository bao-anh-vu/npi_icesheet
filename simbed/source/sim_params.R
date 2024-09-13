## Simulate bed and friction parameters

sim_params <- function(nsims, domain, bed_obs) {
    
    # years <- 20
    # ssa_steady <- steady_state
    # domain <- ssa_steady$domain

    secpera <- 31556926
    fric_scale <- 1e6 * secpera^(1 / 3)

    # 1. Simulate friction coefficients
    print("Simulating friction coefficient...")

    fric.sill <- 8e-5
    fric.nugget <- 0
    fric.range <- 10e3

    simulated_friction <- simulate_friction2(
        nsim = nsims, domain = domain,
        sill = fric.sill, nugget = fric.nugget,
        range = fric.range
    ) 

    simulated_friction <- simulated_friction / fric_scale
    # sim_fric_list <- lapply(1:N, function(c) simulated_friction[, c])

    ## Simulate beds
    print("Simulating beds...")
    
    bed_sims <- simulate_bed(
        nsim = nsims, domain = domain,
        obs_location = bed_obs$locations, obs = bed_obs$obs
    )

    # sim_param_list <- lapply(1:nsims, function(c) {
    #     list(
    #         friction = simulated_friction[, c],
    #         bedrock = bed_sims[, c]
    #     )
    # })

    sim_param_list <- list(
        friction = t(simulated_friction),
        bedrock = t(bed_sims)
    )


    ## DELETE LATER: (just creating an artificial error for debugging purposes)
    # sim_param_list[[2]] <- 0
    # sim_param_list[[3]] <- 0
    
    return(sim_param_list)

}

exp_cov <- function(d, l) {
    return(exp(-3*d / l))
}
