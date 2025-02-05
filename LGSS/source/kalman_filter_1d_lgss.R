kalman_filter_1d_lgss <- function(obs, ini_state_mean, ini_state_var, phi, sigma_eta, sigma_eps) {
    
    # Set up empty vectors to store analysis mean and variance
    x_mean <- c()
    x_var <- c()
    x_mean <- ini_state_mean
    x_var <- ini_state_var

    I <- diag(1, 2, 2)
    H <- 1

    for (k in 2:iters) {
        # Forecast
        x_mean[k] <- phi * x_mean[k-1] #+ rnorm(1, 0, sqrt(tau2))
        P <- phi * x_var[k-1] * phi + sigma_eta^2
        
        # Likelihood as derived by integrating out the states x_t
        mu_x <- P/(sigma_eps^2 + P) * y[k] + sigma_eps^2/(sigma_eps^2 + P) * x_mean[k]
        sigma_x2 <- 1/(1/sigma_eps^2 + 1/P)
        pred_likelihood <- log(1/sqrt(2*pi*(sigma_eps^2 + P)) *
                    exp(1/2 * (mu_x^2 / sigma_x2 - y[k]^2 / sigma_eps^2 - x_mean[k]^2 / P)))
        
        log_likelihood <- log_likelihood + pred_likelihood
        
        # Analysis
        K <- P * t(H) / (H * P * t(H) + sigma_eps^2)
        pseudo_obs <- H * x_mean[k]
        x_mean[k] <- x_mean[k] + K * (y[k] - pseudo_obs)
        x_var[k] <- (1 - K * H) * P
    }

    return(list(x_mean = x_mean, x_var = x_var, log_likelihood = log_likelihood))
}


