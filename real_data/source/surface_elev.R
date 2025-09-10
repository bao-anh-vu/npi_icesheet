get_surface_elev <- function(H, b, z0 = 0, rho = 917.0, rho_w = 1028.0, include_GL = TRUE) {
  GL <- gl_migrate(H, b, z0, rho, rho_w)
  
  z <- c()
  if (include_GL) {
    z[1:GL] <- H[1:GL] + b[1:GL]
    z[(GL+1):length(H)] <- H[(GL+1):length(H)] * (1 - rho/rho_w) + z0
  } else {
    z <- H + b
  }

  z
  
}


gl_migrate <- function (H, b, z0 = 0, rho = 917.0, rho_w = 1028.0) {
  # To calculate GL position given the ice thickness H and bed b
  # Flotation criterion: rho * H >= rho_w * b for grounded ice

  grounded <- rho * H > rho_w * (z0 - b)
  # GL <- sum(grounded)

# plot(rho * H, type = "l")
# lines(rho_w * (z0 - b), col = "red")

  first_false <- which(!grounded)[1]   # index of first FALSE
  if (is.na(first_false)) {
    GL <- length(grounded)                   # if there’s no FALSE, count all
  } else {
    GL <- first_false - 1             # number of TRUEs before it
  }

  return(GL)
}

gl_migrate2 <- function (zs, b, z0 = 0, rho = 917.0, rho_w = 1028.0) {
  # To calculate GL position given the surface elevation zs and bed b
  # Alternative flotation criterion: zs >= (z0 - b) * (rho_w / rho - 1) for grounded ice
  grounded <- zs >= (z0 - b) * (rho_w / rho - 1)
  GL <- sum(grounded, na.rm = TRUE)
  GL
}

calculate_thickness <- function(z, b, z0 = 0, rho = 917.0, rho_w = 1028.0) {
  # 0. Flotation criterion: rho * H >= rho_w * b for grounded ice
  GL <- gl_migrate2(zs = z, b = b)
  H <- z - b
  float_ind <- GL:length(z)
  H[float_ind] <- (z[float_ind] - z0) / (1 - rho / rho_w)
  return(H)
}