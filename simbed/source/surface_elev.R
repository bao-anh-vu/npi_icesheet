surface_elev <- function(H, b, z0 = 0, rho = 910.0, rho_w = 1028.0) {
  include_GL <- TRUE
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


gl_migrate <- function (H, b, z0 = 0, rho = 910.0, rho_w = 1028.0) {
  # To calculate GL position given the ice thickness H and bed b
  # Flotation criterion: rho * H >= rho_w * b for grounded ice
  grounded <- rho * H > rho_w * (z0 - b)
  GL <- sum(grounded)
  
  return(GL)
}

gl_migrate2 <- function (zs, b, z0 = 0, rho = 910.0, rho_w = 1028.0) {
  # To calculate GL position given the surface elevation zs and bed b
  # Alternative flotation criterion: zs >= (z0 - b) * (rho_w / rho - 1) for grounded ice
  grounded <- zs >= (z0 - b) * (rho_w / rho - 1)
  GL <- sum(grounded)
  GL
}

calculate_thickness <- function(z, b, z0 = 0, rho = 910.0, rho_w = 1028.0) {
  # 0. Flotation criterion: rho * H >= rho_w * b for grounded ice
  grounded <- rho * H > rho_w * (z0 - b)
  GL <- sum(grounded)
  GL
}