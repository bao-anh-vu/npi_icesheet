find_gl <- function (thickness, bed, z0 = 0, rho_i = 917.0, rho_w = 1028.0) {
  # To calculate GL position given the ice thickness H and bed b
  # Flotation criterion: rho * H >= rho_w * b for grounded ice, otherwise floating
    H <- thickness
    b <- bed
    grounded <- rho_i * H >= rho_w * (z0 - b)
    GL <- sum(grounded)

    return(GL)
}