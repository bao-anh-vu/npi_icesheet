sample_from_posterior <- function(n, mean, prec_chol) {
    L <- prec_chol #Lmats[[s]]

    ### Sample from the posterior distribution
    z <- matrix(rnorm(S * n_mean_elements), 
                nrow = n_mean_elements, ncol = S)
    v <- backsolve(L, z)

    meanmat <- matrix(rep(mean, S), 
                nrow = n_mean_elements, ncol = S)
    pred_samples <- meanmat + v

    # ## Un-standardise output
    # ## Then dissect the pred samples to get the friction, bed, gl
    # fric_coef_samples <- pred_samples[1:n_fric_basis, ]
    # bed_coef_samples <- pred_samples[(n_fric_basis+1):(n_fric_basis+n_bed_basis), ]
    # gl_samples <- pred_samples[(n_fric_basis+n_bed_basis+1):n_mean_elements, ]

    # fric_coef_samples_ustd <- fric_coef_samples * test_data$sd_fric_coefs + test_data$mean_fric_coefs
    # bed_coef_samples_ustd <- bed_coef_samples * test_data$sd_bed_coefs + test_data$mean_bed_coefs
    # gl_samples_ustd <- gl_samples * test_data$sd_gl + test_data$mean_gl

    # ## Compute friction samples  
    # if (log_transform) {
    #     fric_samples_ls[[s]] <- exp(fric_basis_mat %*% fric_coef_samples_ustd)
    # } else {
    #     fric_samples_ls[[s]] <- fric_basis_mat %*% fric_coef_samples_ustd
    # }

    # bed_samples_ls[[s]] <- bed_basis_mat %*% bed_coef_samples_ustd + bed_mean_mat
    # gl_samples_ls[[s]] <- gl_samples_ustd

    # fric_lq[[s]] <- apply(fric_samples_ls[[s]], 1, quantile, probs = 0.05)
    # fric_uq[[s]] <- apply(fric_samples_ls[[s]], 1, quantile, probs = 0.95)
    
    # bed_lq[[s]] <- apply(bed_samples_ls[[s]], 1, quantile, probs = 0.05)
    # bed_uq[[s]] <- apply(bed_samples_ls[[s]], 1, quantile, probs = 0.95)

    # gl_lq[[s]] <- apply(gl_samples_ls[[s]], 1, quantile, probs = 0.05)
    # gl_uq[[s]] <- apply(gl_samples_ls[[s]], 1, quantile, probs = 0.95)
}