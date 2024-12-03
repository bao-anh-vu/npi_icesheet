cond_sim_gp <- function(nsims, x_test, x_train, obs, obs_sd, mu_obs) {

    eps <- 1e-6
    # x_test <- matrix(seq(0, 2*pi, length=n), ncol=1)
    y <- obs

    mu <- get_bed_mean(obs_locations = x_train, obs = obs, 
                        domain = sort(c(x_train, x_test)))
                    
    # mu_y <- mu %>% filter(x %in% x_train) %>% select(mean) %>% unlist()
    mu_y <- y #mu$mean[mu$x %in% x_train]
    mu_test <- mu$mean[mu$x %in% x_test]

    # if (use_bed_example) {
        Sigma <- get_bed_cov(si = x_train, sj = x_train)#
        S_test_test <- get_bed_cov(si = x_test, sj = x_test)#
        S_test_train <- get_bed_cov(si = x_test, sj = x_train)
    # } else {
        # D <- rdist(x_train) 
        # Dx_test <- rdist(x_test)
        # Dx_test_train <- rdist(x_test, x_train)
        
    #     Sigma <- exp(-D) + diag(eps, ncol(D))
    #     S_test_test <- exp(-Dx_test) + diag(eps, ncol(Dx_test))
    #     S_test_train <- exp(-Dx_test_train) 
    # }

    Si <- solve(Sigma + diag(obs_sd, ncol(Sigma))) 
    mu_p <- mu_test + S_test_train %*% Si %*% (y - mu_y)
    Sigma_p <- S_test_test - S_test_train %*% Si %*% t(S_test_train)

    sims <- rmvnorm(nsims, mu_p, Sigma_p)

    # if (!is.null(obs_sd)) {
    #     sims_at_train <- t(rmvnorm(nsims, mu_y, diag(obs_sd^2)))
    # } else {
    #     sims_at_train <- matrix(rep(mu_y, nsims), ncol=nsims)
    # }

    # test_df <- cbind(x_test, mu_p, t(sims))
    # train_df <- cbind(x_train, mu_y, sims_at_train)
    # sim_df <- as.data.frame(rbind(test_df, train_df))
    # colnames(sim_df) <- c("x", "mean", paste0("sim", 1:nsims))

    # sim_df <- dplyr::arrange(sim_df, x)
    sim_df <- data.frame(x_test, sims)
    colnames(sim_df) <- c("x", paste0("sim", 1:nsims))

    return(list(mu_p = mu_p, Sigma_p = Sigma_p, sims = sim_df))
}

get_bed_cov <- function(si, sj, sill = 4000, nugget = 200, range = 50e3) { 
  # matrix of distances
  d <- rdist(si, sj)

  # nugget model
  nug <- ifelse(d > 0, nugget, 0)

  # variogram
  psill <- sill - nugget
  variogram <- psill * (1 - exp(-3 * d / range)) + nug

  #covariance
  cov <- sill - variogram
  # s * (1 - exp(-3 * spDists(coordinates(x),coordinates(y)) / range))
  return(cov)
}

get_bed_mean <- function(obs_locations, obs, domain) {
    ## Need to account for the fact that the bed has a non-zero mean
    df <- data.frame(obs_locations = obs_locations, bed_elev = obs)
  
    bed_fit <- loess(bed_elev ~ obs_locations, data = df, span = 0.4,
                    control = loess.control(surface = "direct")) 
    # bed.fit <- lm(bed_elev ~ poly(obs_locations, 10), data = df) #, span = 0.12,

    ## Use loess to interpolate between observations
    ## Outside of that range just use the first and last values in the fit
    lower_tail <- which(domain < obs_locations[1])
    upper_tail <- which(domain > obs_locations[length(obs_locations)])
    bed_mean <- predict(bed_fit, newdata = data.frame(obs_locations = domain))
    bed_mean[lower_tail] <- bed_fit$fitted[1]
    bed_mean[upper_tail] <- bed_fit$fitted[length(bed_fit$fitted)]

    mean_df <- data.frame(x = domain, mean = bed_mean)

    return(mean_df)
}
