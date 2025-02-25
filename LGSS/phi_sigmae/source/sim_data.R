## LGSS forward model

sim_data <- function(phi, sigma_eta = 0.7, sigma_eps = 0.5, iters) {
    
    ## Simulate the states
    x <- c()
    x[1] <- rnorm(1, 0, sqrt(sigma_eta^2 / (1 - phi^2)))
    
    for (t in 2:iters) {
        x[t] <- phi * x[t-1] + rnorm(1, 0, sigma_eta)
    }

    # Generate observations y_1:T
    y <- x + rnorm(iters, 0, sigma_eps)

    return(data = list(x = x, y = y))
}