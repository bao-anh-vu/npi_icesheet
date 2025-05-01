solve_velocity <- function(p, L, J, H, b, C, u, u_bc) {

  # % SSAFLOWLINE  Computes the velocity from the SSA in a flow-line case, with
  # %   left-end Dirichlet condition  u=ug  and right-end Neumann condition at
  # %   calving-front.
  # % form:
  # %   [u,u0] = ssaflowline(p,J,H,b,uleft,initchoice)
  # % all 6 input arguments are required:
  # %   p  = list with parameters A,C,g,L,m,n,rho,rhow,secpera
  # %   J  = number of grid subintervals
  # %   H  = thickness (m), a length J+1 column vector
  # %   b  = bed elevation (m), same size as H
  # %   u_bc = velocity boundary value at x=0 end of domain
  # %   initchoice = 1,2 chooses initial guess scheme;
  # %                use 1 for ice shelves and 2 for ice streams; NOT USED ATM
  # % outputs:
  # %   u = the velocity solution (m s-1), a length J+1 column vector
  # %   u0 = the velocity initial guess (m s-1), same size as u
  # % Does Picard iteration to solve nonlinear SSA problem.
  # % Calls:  FLOWLINE to solve inner linear PDE boundary value problem,
  # %         SSAINIT  to get initial iterate

  # Check if all required arguments are provided
  if (missing(p) || missing(J) || missing(H) || missing(u) || missing(b) || missing(C) || missing(u_bc)) {
    stop("All 6 input arguments are required.")
  }

  dx <- L / J
  x <- seq(0, L, length.out = J + 1)
  xstag <- seq(dx/2, L + dx/2, by = dx)

  ## Determine GL position
  GL <- find_gl(H, b, rho_i = p$rho_i, rho_w = p$rho_w)
  

  # Drag coefficient
  alpha <- rep(0, length(x)) # initialize alpha
  alpha[1:GL] <- C[1:GL] * abs(u[1:GL])^(p$m-1) # for grounded ice
  # alpha <- C

  # Ice rigidity

  # Surface elevation (h here is z_s in my notation)
  zs <- H + b # for grounded ice
  zs[(GL+1):(J+1)] <- (1 - p$rho_i/p$rho_w) * H[(GL+1):(J+1)] # for floating ice
  
  # Compute regular grid values of slope (dz_s/dx)
  dzdx <- numeric(length(H))
  dzdx[1] <- (zs[2] - zs[1]) / dx
  dzdx[2:J] <- (zs[3:(J+1)] - zs[1:(J-1)]) / (2 * dx)
  dzdx[J+1] <- (zs[J+1] - zs[J]) / dx

  # Right hand side of linearised system
  beta <- p$rho_i * p$g * H * dzdx

  # Boundary condition: velocity at the calving front
  gamma <- (0.25 * p$A^(1/p$n) * (1 - p$rho_i/p$rho_w) * p$rho_i * p$g * H[J+1])^p$n
  # gamma <- (0.25 * p$B^(-1) * (1 - p$rho_i/p$rho_w) * p$rho_i * p$g * H[J+1])^p$n

  # Initial velocity (modify later to use ssainit())
  # u0 <- 0.001 * x
  # u0[1] <- u_bc # boundary condition at the ice divide (x1 = 0)
  # # u0 <- ssainit(p, x, beta, gamma, initchoice)
  u0 <- u

  Hstag <- 0.5 * (H[1:J] + H[2:(J+1)])
  tol <- 1 / p$secpera #1.0e-14
  eps_reg <- (1.0 / p$secpera) / L
  maxdiff <- Inf
  W <- numeric(J + 1)
  iter <- 0

  while (maxdiff > tol) {
    uxstag <- (u[2:(J+1)] - u[1:J]) / dx
    sqr_ux_reg <- uxstag^2 + eps_reg^2

    W[1:J] <- 2 * p$A^(-1/p$n) * Hstag * sqr_ux_reg^(((1/p$n) - 1) / 2)
    # W[1:J] <- 2 * p$B * Hstag * sqr_ux_reg^(((1/p$n) - 1) / 2)
    W[J+1] <- W[J]

    unew <- flowline(L, J, gamma, W, alpha, beta, u_bc, GL)
    maxdiff <- max(abs(unew - u))
    u <- unew
    iter <- iter + 1
    # cat(".")

  }

  # cat(sprintf("\nSSA solver did %d Picard iterations on dx = %.3f km grid\n", iter, dx / 1000))
  return(list(u = u, u0 = u0, zs = zs, b = b, H = H, GL = GL))

}
