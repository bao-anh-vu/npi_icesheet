flowline_eqns <- function(huxg, bed, params) {
  Nx <- params$Nx
  N1 <- params$N1
  ds <- params$dsigma
  sigma <- params$sigma
  sigma_elem <- params$sigma_elem
  
  # Unpack h, u, xg
  h <- huxg[1:Nx]
  u <- huxg[(Nx + 1):(2 * Nx)]
  xg <- huxg[2 * Nx + 1]
  
  if (use_linear_bed) {
    # Floatation thickness
      hf <- (-bed(xg * params$xscale, params) / params$hscale) / (1 - params$lambda)
    
    # Bed profile
    b <- -bed(xg * sigma * params$xscale, params) / params$hscale

    ## Not sure why there is a minus sign in front of the bed for Alex's solver but that's how it works

  } else {
    b <- bed / params$hscale
    hf <- b[which(x == xg)] / (1 - params$lambda)
    
  }

  # Derived values
  m <- params$m
  nglen <- params$nglen
  lambda <- params$lambda
  a <- params$accum / params$ascale
  eps <- params$eps
  
  # Mass conservation (Fh)
  Fh <- numeric(Nx)
  Fh[1] <- (2 * h[1] * u[1]) / (ds[1] * xg) - a
  Fh[2] <- (h[2] * (u[2] + u[1])) / (2 * xg * ds[2]) - a
  for (i in 3:(Nx - 1)) {
    term1 <- h[i] * (u[i] + u[i - 1])
    term2 <- h[i - 1] * (u[i - 1] + u[i - 2])
    Fh[i] <- (term1 - term2) / (2 * xg * ds[i]) - a
  }
  Fh[N1] <- (1 + 0.5 * (1 + ds[N1] / ds[N1 - 1])) * h[N1] -
            0.5 * (1 + ds[N1] / ds[N1 - 1]) * h[N1 - 1] -
            h[N1 + 1]
  Fh[Nx] <- (h[Nx] * (u[Nx] + u[Nx - 1]) -
             h[Nx - 1] * (u[Nx - 1] + u[Nx - 2])) / (2 * xg * ds[Nx - 1]) - a
  
  # Momentum conservation (Fu)
  Fu <- numeric(Nx)
  Fu[1] <- (4 * eps) * (1 / (xg * ds[1])^(1 / nglen + 1)) * (
    h[2] * (u[2] - u[1]) * abs(u[2] - u[1])^(1 / nglen - 1) -
    h[1] * (2 * u[1]) * abs(2 * u[1])^(1 / nglen - 1)
  ) -
    u[1] * abs(u[1])^(m - 1) -
    0.5 * (h[1] + h[2]) * (h[2] - b[2] - h[1] + b[1]) / (xg * ds[1])
  
  for (i in 2:(Nx - 1)) {
    Fu[i] <- (4 * eps) * (1 / (xg * ds[i])^(1 / nglen + 1)) * (
      h[i + 1] * (u[i + 1] - u[i]) * abs(u[i + 1] - u[i])^(1 / nglen - 1) -
      h[i] * (u[i] - u[i - 1]) * abs(u[i] - u[i - 1])^(1 / nglen - 1)
    ) -
      u[i] * abs(u[i])^(m - 1) -
      0.5 * (h[i] + h[i + 1]) * (h[i + 1] - b[i + 1] - h[i] + b[i]) / (xg * ds[i])
  }
  
  Fu[N1] <- (u[N1 + 1] - u[N1]) / ds[N1] - (u[N1] - u[N1 - 1]) / ds[N1 - 1]

  Fu[Nx] <- (1 / (xg * ds[Nx - 1])^(1 / nglen)) *
    (abs(u[Nx] - u[Nx - 1])^(1 / nglen - 1)) * (u[Nx] - u[Nx - 1]) -
    lambda * hf / (8 * eps)
  
  # GL condition
  Fxg <- 3 * h[Nx] - h[Nx - 1] - 2 * hf
  
  return(c(Fh, Fu, Fxg))
}
