cond_sim <- function(mu, sd, taxis, l, nsims = 1L) {
  
  stopifnot(length(mu) == length(sd))
  stopifnot(length(mu) == length(taxis))
  stopifnot(length(l) == 1L & length(nsims) == 1L) 
  N <- length(mu)
  
  ## Create GP assuming e-folding length scale of 5 years
  D <- as.matrix(dist(taxis))
  # C <- outer(sd, sd) * exp(-D/l)
  C <- outer(sd, sd) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
  L <- t(chol(C))
  simfields <- matrix(mu, ncol = nsims, nrow = length(mu), byrow = FALSE) + 
    L %*% matrix(rnorm(N*nsims), ncol = nsims) %>% as.data.frame()
  names(simfields) <- paste0("Sim", 1:nsims)
  simfields$taxis <- taxis
  simfields
}

cond_sim_sigma <- function(sd, taxis, l) {
  ## Create GP assuming e-folding length scale of 5 years
  D <- as.matrix(dist(taxis))
  # C <- outer(sd, sd) * exp(-D/l)
  C <- outer(sd, sd) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
  return(C)
}