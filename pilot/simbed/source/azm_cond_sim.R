cond_sim <- function(mu, sd, taxis, l, nsims = 1L) {
  
  stopifnot(length(mu) == length(sd))
  stopifnot(length(mu) == length(taxis))
  stopifnot(length(l) == 1L & length(nsims) == 1L) 
  N <- length(mu)
  
  ## Create GP assuming e-folding length scale of 5 years
  # D <- as.matrix(dist(taxis))
  D <- rdist(taxis)
  
  C <- outer(sd, sd) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
  L <- t(chol(C))
  
  # simfields <- matrix(mu, ncol = nsims, nrow = length(mu), byrow = FALSE) + 
  #   L %*% matrix(rnorm(N*nsims), ncol = nsims) %>% as.data.frame()
  simfields <- matrix(mu, ncol = nsims, nrow = length(mu), byrow = FALSE) +
    L %*% matrix(rnorm(N*nsims), ncol = nsims)

  # simfields <- as.data.frame(simfields_mat)
  # names(simfields) <- paste0("Sim", 1:nsims)
  # simfields$taxis <- taxis

  ## Generate noise using the correlation matrix instead of covariance
  # ones <- rep(1, length(taxis))
  # R <- outer(ones, ones) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
  # L2 <- t(chol(R))
  # C2 <- diag(sd) %*% R %*% diag(sd)

  
  # # then later do
  # sim <- mu + diag(sd) %*% L2 %*% rnorm(N, 0, 1)

  simfields
}

cond_sim2 <- function(mu, sd, L, nsims = 1L) {
  
  stopifnot(length(mu) == length(sd))
  # stopifnot(length(mu) == length(taxis))
  # stopifnot(length(l) == 1L & length(nsims) == 1L) 
  N <- length(mu)
  sd <- as.vector(sd)
  ## Create GP assuming e-folding length scale of 5 years
  # L <- Matrix(L, sparse = T)
  # L <- as(L, "dgCMatrix")
  simfields <- mu + sd * (L %*% rnorm(N, 0, 1))
  
  # simfields <- matrix(mu, ncol = nsims, nrow = length(mu), byrow = FALSE) +
    # L %*% matrix(rnorm(N*nsims), ncol = nsims)
  simfields
}

cond_sim_sigma <- function(sd, taxis, l) {
  ## Create GP assuming e-folding length scale of 5 years
  D <- as.matrix(dist(taxis))
  # C <- outer(sd, sd) * exp(-D/l)
  C <- outer(sd, sd) * (1 + sqrt(3) * D / l) * exp(- sqrt(3) * D / l )
  return(C)
}