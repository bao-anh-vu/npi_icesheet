process_sim_results <- function(sims, sim_beds = T) {

    N <- length(sims)
    
    ## Note to self: should save a version with the "true" friction coef as well
    # thickness_velocity_list <- lapply(sims, function(sim) sim$thickness_velocity_arr)
    surface_obs_list <- lapply(sims, function(sim) sim$surface_obs)
    friction_list <- lapply(sims, function(sim) sim$friction_arr)
    bed_list <- lapply(sims, function(sim) sim$bed_arr)
    gl_list <- lapply(sims, function(sim) sim$grounding_line)

    true_surface_elevs <- lapply(sims, function(sim) array(sim$true_surface_elevs, dim = c(dim(sim$true_surface_elevs), 1)))
    true_thicknesses <- lapply(sims, function(sim) array(sim$true_thicknesses, dim = c(dim(sim$true_thicknesses), 1)))
    true_velocities <- lapply(sims, function(sim) array(sim$true_velocities, dim = c(dim(sim$true_velocities), 1)))

    # concat_input <- do.call("rbind", thickness_velocity_list)
    # thickness_velocity_arr <- array(concat_input, dim = c(N, dim(thickness_velocity_list[[1]])))
    concat_input <- do.call("rbind", surface_obs_list)
    surface_obs_arr <- array(concat_input, dim = c(N, dim(surface_obs_list[[1]])))

    concat_output <- do.call("rbind", friction_list)
    friction_arr <- array(concat_output, dim = c(N, length(friction_list[[1]])))

    if (sim_beds) {
        concat_output <- do.call("rbind", bed_list)
        bed_arr <- array(concat_output, dim = c(N, length(bed_list[[1]])))
    }

    concat_gl <- do.call("rbind", gl_list)
    gl_arr <- array(concat_gl, dim = c(N, length(gl_list[[1]])))

    concat_true_se <- do.call("rbind", true_surface_elevs)
    true_se_arr <- array(concat_true_se, dim = c(N, dim(true_surface_elevs[[1]])))

    concat_true_th <- do.call("rbind", true_thicknesses)
    true_th_arr <- array(concat_true_th, dim = c(N, dim(true_thicknesses[[1]])))

    concat_true_v <- do.call("rbind", true_velocities)
    true_v_arr <- array(concat_true_v, dim = c(N, dim(true_velocities[[1]])))

    # if (sim_beds) {
        return_obj <- list(
            # thickness_velocity_arr = thickness_velocity_arr,
            surface_obs_arr = surface_obs_arr,
            friction_arr = friction_arr,
            bed_arr = bed_arr,
            gl_arr = gl_arr,
            true_surface_elevs = true_se_arr,
            true_thicknesses = true_th_arr,
            true_velocities = true_v_arr
        )
    # } else {
    #     return_obj <- list(
    #         # thickness_velocity_arr = thickness_velocity_arr,
    #         surface_obs_arr = surface_obs_arr,
    #         friction_arr = friction_arr,
    #         gl_arr = gl_arr
    #     )
    # }
    return(return_obj)
}