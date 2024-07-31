fit_friction_basis <- function(nbasis, domain, sample_arr) {
    ## Place basis function centres along domain
    # basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis+2)
    # basis_centres <- basis_centres[2:(length(basis_centres)-1)] 
                            
    basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis)
    testbasis <- local_basis(manifold = real_line(), 
                            loc = matrix(basis_centres),
                            type = "bisquare",
                            scale = rep(5e3, nbasis))
    # show_basis(testbasis)
    basis_fns <- lapply(testbasis@fn, function(f) f(domain))
    # plot(domain, basis_fns[[10]])
    basis_mat <- do.call(cbind, basis_fns)
    # matplot(domain/1000, basis_mat, type = "l", col= "salmon", lty = 1, lwd = 1.5,
    #         xlab = "Domain (km)", xlim = c(0, 200))


    N <- dim(sample_arr)[1]
    basis_coefs <- matrix(NA, N, nbasis) 
    fitted_values <- matrix(NA, N, length(domain))
    for (sim in 1:N) { # parallelise
        df_local <- as.data.frame(cbind(sample_arr[sim, ], basis_mat))
        colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
        lmfit_local <- lm(fric ~ . - 1, data = df_local) # no intercept for now so -1
        basis_coefs[sim, ] <- as.vector(lmfit_local$coefficients)
        fitted_values[sim, ] <- as.vector(lmfit_local$fitted.values)
    }

    return(list(basis_coefs = basis_coefs, 
                basis_mat = basis_mat,
                fitted_values = fitted_values))

}


## Basis for the bedrock
fit_bed_basis <- function(nbasis, domain, sample_arr) {
    ## Place basis function centres along domain
    # basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis+2)
    # basis_centres <- basis_centres[2:(length(basis_centres)-1)] 
                            
    basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis)
    testbasis <- local_basis(manifold = real_line(), 
                            loc = matrix(basis_centres),
                            type = "bisquare",
                            scale = rep(5e3, nbasis))
    # show_basis(testbasis)
    basis_fns <- lapply(testbasis@fn, function(f) f(domain))
    # plot(domain, basis_fns[[10]])
    basis_mat <- do.call(cbind, basis_fns)
    # matplot(domain/1000, basis_mat, type = "l", col= "salmon", lty = 1, lwd = 1.5,
    #         xlab = "Domain (km)", xlim = c(0, 200))

    N <- dim(sample_arr)[1]

    # test1 <- system.time({
    #     basis_coefs <- matrix(NA, N, nbasis)
    #     fitted_values <- matrix(NA, N, length(domain))
    #     for (sim in 1:N) { # parallelise
    #         df_local <- as.data.frame(cbind(sample_arr[sim, ], basis_mat))
    #         colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
    #         lmfit_local <- lm(fric ~ . - 1, data = df_local) # no intercept for now so -1
    #         basis_coefs[sim, ] <- as.vector(lmfit_local$coefficients)
    #         fitted_values[sim, ] <- as.vector(lmfit_local$fitted.values)
    #     }
    # }
    # )

    ## Parallel ver
    # test2 <- system.time({
        basis_fit <- mclapply(1:N, function(sim) {
        df_local <- as.data.frame(cbind(sample_arr[sim, ], basis_mat))
        colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
        lmfit_local <- lm(fric ~ . - 1, data = df_local) # no intercept for now so -1
        coefs <- as.vector(lmfit_local$coefficients)
        fitted_values <- as.vector(lmfit_local$fitted.values)
        return(list(coefs = coefs,
                    fitted_values = fitted_values))
        },
        mc.cores = 20L
    )
    basis_coefs <- do.call(rbind, lapply(basis_fit, function(x) x$coefs))
    fitted_values <- do.call(rbind, lapply(basis_fit, function(x) x$fitted_values))
    # })
   
# plot(sample_arr[1,], type = "l")
# lines(fitted_values[1, ], col = "red")
#     browser()

    return(list(basis_coefs = basis_coefs, 
                basis_mat = basis_mat,
                fitted_values = fitted_values))

}