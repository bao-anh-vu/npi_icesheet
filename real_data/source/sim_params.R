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
    
    bed_sim_output <- simulate_bed(N, domain = domain, 
                            obs_locations = bed_obs$ind, 
                            obs = bed_obs$bed_elev, 
                            sill = 1000, nugget = 0) 

    bed_sims <- bed_sim_output$sims
    bed_mean <- bed_sim_output$mean

    # out <- cond_sim_gp(1, x_test = bed_obs$loc[bed_obs$chosen == 0], 
    #                     x_train = bed_obs$loc[bed_obs$chosen == 1],
    #                     obs = bed_obs$bed_elev[bed_obs$chosen == 1], 
    #                     obs_sd = bed_obs$bed_sd[bed_obs$chosen == 1])
    # bed_sims <- out$sims[-(1:2), ]

    sim_param_list <- list(
        friction = t(simulated_friction),
        bedrock = t(bed_sims)
    )

    ## DELETE LATER: (just creating an artificial error for debugging purposes)
    # sim_param_list[[2]] <- 0
    # sim_param_list[[3]] <- 0
    
    return(list(sim_param_list = sim_param_list, bed_mean = bed_mean))

}

exp_cov <- function(d, l) {
    return(exp(-3*d / l))
}
