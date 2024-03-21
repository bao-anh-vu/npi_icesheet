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
    basismat <- do.call(cbind, basis_fns)
    # matplot(domain/1000, basismat, type = "l", col= "salmon", lty = 1, lwd = 1.5,
    #         xlab = "Domain (km)", xlim = c(0, 200))


    N <- ncol(friction_arr)
    fitted_fric <- matrix(NA, length(domain), N)
    for (sim in 1:N) {
    df_local <- as.data.frame(cbind(friction_arr[, sim], basismat))
    colnames(df_local) <- c("fric", sapply(1:nbasis, function(x) paste0("eof", x)))
    lmfit_local <- lm(fric ~ ., data = df_local)
    fitted_fric[, sim] <- lmfit_local$fitted.values
    }

    if (save_sims) {
        saveRDS(fitted_fric, file = paste0("./output/fitted_fric_", simset, "_", data_date))
    }

    return(fitted_fric)

}