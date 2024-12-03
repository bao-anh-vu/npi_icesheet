get_ini_thickness <- function(surf_elev, bed, rho = 910, rho_w = 1028) {
    se_grounded <- surf_elev # data is only available where ice is grounded 
    H_ini_ground <- se_grounded - bed[1:length(se_grounded)]
    missing <- which(is.na(surf_elev)) # missing data is where ice is floating
    # rho <- 910.0
    # rho_w <- 1028.0
    thickness_at_tail <- - bed[missing] * rho_w / rho # work out ice thickness for the floating shelf based on flotation condition
    H_ini_all <- c(na.omit(H_ini_ground), thickness_at_tail) #+ offset
    
    # H_ini_new <- H_ini_new + offset
    return(H_ini_all)
}