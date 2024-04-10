fit_fric_basis <- function(nbasis, domain, friction_arr) {
    basis_centres <- seq(domain[1], domain[length(domain)], length.out = nbasis+2)
    basis_centres <- basis_centres[2:(length(basis_centres)-1)]
                            
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


    N <- dim(friction_arr)[1]
    basis_coefs <- matrix(NA, N, nbasis+1) # +1 for the intercept
    fitted_fric <- matrix(NA, N, length(domain))
    for (sim in 1:N) {
        df_local <- as.data.frame(cbind(friction_arr[sim,,,], basis_mat))
        colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
        lmfit_local <- lm(fric ~ ., data = df_local)
        basis_coefs[sim, ] <- as.vector(lmfit_local$coefficients)
        fitted_fric[sim, ] <- as.vector(lmfit_local$fitted.values)
    }

    return(list(basis_coefs = basis_coefs, 
                basis_mat = basis_mat,
                fitted_fric = fitted_fric))

}