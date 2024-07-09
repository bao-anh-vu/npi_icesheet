process_sim_results <- function(sims, sim_beds = T) {

    N <- length(sims)
    
    ## Note to self: should save a version with the "true" friction coef as well
    thickness_velocity_list <- lapply(sims, function(sim) sim$thickness_velocity_arr)
    friction_list <- lapply(sims, function(sim) sim$friction_arr)
    bed_list <- lapply(sims, function(sim) sim$bed_arr)
    gl_list <- lapply(sims, function(sim) sim$grounding_line)

    concat_input <- do.call("rbind", thickness_velocity_list)
    thickness_velocity_arr <- array(concat_input, dim = c(N, dim(thickness_velocity_list[[1]])))

    concat_output <- do.call("rbind", friction_list)
    friction_arr <- array(concat_output, dim = c(N, length(friction_list[[1]])))

    if (sim_beds) {
        concat_output <- do.call("rbind", bed_list)
        bed_arr <- array(concat_output, dim = c(N, length(bed_list[[1]])))
    }

    concat_gl <- do.call("rbind", gl_list)
    gl_arr <- array(concat_gl, dim = c(N, length(gl_list[[1]])))

    if (sim_beds) {
        return_obj <- list(
            thickness_velocity_arr = thickness_velocity_arr,
            friction_arr = friction_arr,
            bed_arr = bed_arr,
            gl_arr = gl_arr
        )
    } else {
        return_obj <- list(
            thickness_velocity_arr = thickness_velocity_arr,
            friction_arr = friction_arr,
            gl_arr = gl_arr
        )
    
    }
    return(return_obj)
}