solve_velocity <- function(prev_velocity, thickness, domain, bed, friction, 
                           perturb_hardness = FALSE, include_GL = TRUE, 
                           B_variable = FALSE, velocity_bc = 0, fixed_u0 = TRUE, 
                           secpera = 31556926,
                           n = 3, m = 1/3, rho = 910.0, rho_w = 1028.0,
                           g = 9.81, A = 1.4579e-25, 
                           Bg = 0.4 * 1e6 * 31556926 ^ (1/3), z0 = 0) {
  
  ## Define physical parameters -- pass as "default" function values instead?
  # secpera <- 31556926 #seconds per annum
  # n <- 3.0 # exponent in Glen's flow law
  # m <- 1/n # friction exponent
  # rho <- 910.0 #ice density
  # rho_w <- 1028.0 #sea water density
  # r <- rho / rho_w # define this ratio for convenience
  # g <- 9.81 #gravity constant
  # A <- 1.4579e-25 #ice hardness
  # z0 <- 0 # ocean surface elevation
  u <- prev_velocity
  H <- thickness
   
  u0 <- velocity_bc
  
  ## Convert input velocity into m/s
  u <- u / secpera
  
  ## Domain
  x <- domain
  J <- length(domain)-1 # number of grid 'spaces'
  L <- x[length(x)] # length of flowline
  dx <- x[2] - x[1] # grid spacing
  
  C <- friction 
  b <- bed
  
  ## The equation we're trying to solve is of the form
  ## W(x) du/dx + gamma(x) u = beta(x)
  ## Set up some vectors for gamma, beta, W
  
  gamma <- rep(0, length(x))
  W <- rep(0, length(x))
  beta <- rep(0, length(x))
  
  ## Define a strain rate for regularisation of du/dx
  eps_reg <- (1.0 / secpera) / L # strain rate of 1 m/a over length of shelf
  
  # 1. Compute du_dx on a grid and regularise
  du_dx <- c( t((u[2] - u[1])/dx), t((u[3:(J+1)] - u[1:(J-1)])/(2*dx)),
              t((u[J+1] - u[J])/dx))
  
  du_dx_reg = sqrt(du_dx^2 + eps_reg^2) #regularise to avoid division by zero
  
  # 2. Average values of H on a staggered grid (for stability)
  H_stag <- 0.5 * (H[2:length(H)] + H[1:(length(H)-1)])
  
  if (include_GL) {
    
    GL <- gl_migrate(H, b)

    # cat("GL position: ",  GL / J * L / 1000,  "\n")
    ## Ice hardness
    if (perturb_hardness) {
      Bg <- 0.3 * 1e6 * secpera ^ m
      # Bg <- 0.25 * 1e6 * secpera ^ m
    }
    B <- rep(Bg, length(x))
    #print(Bg)
    
    if (B_variable) {
      # Calculate T(xg)
      T0 <- 0.5 * (1-r) * rho * g * H[GL]^2
      # Calculate B(x) the ice hardness
      B[1:GL] <- T0 / (2 * H[1:GL] * abs(du_dx_reg[1:GL])^(1/n - 1) * du_dx_reg[1:GL])
    }
    
  }
  
  stopifnot(exprs = {
                      !is.null(GL)
                      GL > 1
                    })
  
  if (GL > 1) {
    B_stag <- 0.5 * (B[2:(J+1)] + B[1:J])
    
    # and calculate W(x) using the current iteration of du/dx
    W[1:J] <- 2 * B_stag[1:J] * H_stag * abs(du_dx_reg[1:J])^((1/n)-1)
    W[J+1] <- W[J] # technically these are W_{J+3/2} and W_{J+1/2}
    
    # 3. Construct gamma
    if (include_GL) {
      # gamma[1:GL] <- C[1:GL] * rho * g * H[1:GL]
      gamma[1:GL] <- C[1:GL] * abs(u[1:GL])^(m-1)
      gamma[(GL+1):(J+1)] <- 0
    } else{
      # gamma <- C * rho * g * H
      gamma <- C * abs(u)^(m-1)
    }
    
    # 3. Construct matrix A
    
    ##### CHANGED HERE #######
    ##### CONSTRUCTING A_mat IN ONE GO ########
    ##### COMMENTED OUT ALL OTHER A_mat BELOW #########
    
    #A_mat <- matrix(0, J + 1, J + 1)
    
    v1 <- c(0, W[1 : J-1],W[J] + W[J + 1] ) # first element is dummy
    v2 <- c(0, -(W[1:(J-1)] + W[2:J] + gamma[2:J] * dx^2), -(W[J] + W[J + 1] + dx^2 * gamma[J + 1]))
    v3 <- c(0, W[2 : J], 0) # last element is dummy
    
    if (include_GL) {
      if (GL <= (J - 1)) {
        #A_mat[GL, (GL-1):(GL+1)] <- c(-1, 2, -1)
        v1[GL] <- -1
        v2[GL] <- 2
        v3[GL] <- -1
      }
    }
    
    if (fixed_u0) {
      #A_mat[1, 1] <- 1
      v2[1] <- 1
    } else {
      #A_mat[1, 1:2] <- c(-1, 1)
      v2[1] <- -1
      v3[1] <- 1 
    }
    
    V <- cbind(v1, v2, v3)
    
    row_max <- as.vector(rowMaxs(abs(V))) #find maximum of each row
    V[1:nrow(V),] <- V[1:nrow(V),] / row_max[1:nrow(V)]
    
    A_mat <- bandSparse(J + 1, J + 1, #dimensions
                        (-1):1, #band, diagonal is number 0
                        list(V[2:(J+1),1], V[,2], V[(1:J),3]))
    
    #for (j in 2:J) {
    #  A_mat[j, (j-1):(j+1)] <- c(W[j - 1], -(W[j-1] + W[j] + gamma[j] * dx^2), W[j])
    # }
    
    # 4. Construct vector b on the rhs
    # First compute dz/dx
    z <- c() # top surface elevation
    if (include_GL) {
      z[1:GL] <- H[1:GL] + b[1:GL]
      z[(GL+1):(J+1)] <- H[(GL+1):(J+1)] * (1 - rho / rho_w) + z0
    } else {
      z <- H + b
    }
    
    dz_dx <- c( t((z[2] - z[1])/dx), t((z[3:(J+1)] - z[1:(J-1)])/(2*dx)),
                t((z[J+1] - z[J])/dx))
    
    # then compute the whole rhs
    beta <- as.vector(rho * g * H * dz_dx)
    rhs <- dx^2 * beta
    
    # BC at calving front
    # Last row in A is special
    #A_mat[J + 1, J] <- W[J] + W[J + 1] #this is equal to 2 W[J]
    #A_mat[J + 1, J + 1] <- -(W[J] + W[J + 1] + dx^2 * gamma[J + 1])
    u_c <- ( 0.25 * (1 / Bg) * (1 - rho / rho_w) * rho * g * H[J + 1] )^n
    rhs[J + 1] <- rhs[J + 1] - 2 * u_c * dx * W[J + 1]
    
    # To make the solution continuous,
    # take u(xg + 1) as the average of u(xg + 2) and u(xg)
    if (include_GL) {
      if (GL <= (J - 1)) {
        #A_mat[GL, (GL-1):(GL+1)] <- c(-1, 2, -1)
        rhs[GL] <- 0
      }
    }
    
    # Boundary condition at the ice divide (x = 0)
    if (fixed_u0) {
      #A_mat[1, 1] <- 1
      rhs[1] <- u0
    } else {
      #A_mat[1, 1:2] <- c(-1, 1)
      rhs[1] <- 0
    }
    
    # 5. Scaling to avoid singularity
    # row_max <- apply(abs(A_mat), 1, max) #find maximum of each row
    
    #row_max <- as.vector(rowMax(abs(A_mat))) #find maximum of each row
    #A_mat[1:nrow(A_mat),] <- A_mat[1:nrow(A_mat),] / row_max[1:nrow(A_mat)]
    rhs <- rhs / row_max
    
    #A_sparse <- as(A_mat, "dgCMatrix")
    A_sparse <- A_mat
    
    # 6. Solve for the velocity
    u_new <- solve(A_sparse, rhs)
    
  } else {
    u_new <- u
  }
  
  ## Convert the velocity into m/yr again
  u_new <- u_new * secpera
  
  return(u_new)
  
}
