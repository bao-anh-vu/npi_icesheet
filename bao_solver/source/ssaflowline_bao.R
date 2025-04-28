solve_velocity <- function(p, J, H, b, C, u, u_bc) {

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

  dx <- p$L / J
  x <- seq(0, p$L, length.out = J + 1)
  xstag <- seq(dx/2, p$L + dx/2, by = dx)

  ## Determine GL position
  GL <- find_gl(H, b, rho_i = param$rho, rho_w = param$rhow)

  # Drag coefficient
  # alpha <- 0 # initialize alpha
  # alpha[1:GL] <- rep(p$C, length(x)) # for grounded ice

  alpha <- C

  # Surface elevation (h here is z_s in my notation)
  h <- H + b # for grounded ice
  h[(GL+1):(J+1)] <- (1 - p$rho/p$rhow) * H[(GL+1):(J+1)] # for floating ice
  
  # Compute regular grid values of slope (dz_s/dx)
  hx <- numeric(length(H))
  hx[1] <- (h[2] - h[1]) / dx
  hx[2:J] <- (h[3:(J+1)] - h[1:(J-1)]) / (2 * dx)
  hx[J+1] <- (h[J+1] - h[J]) / dx

  # Right hand side of linearised system
  beta <- p$rho * p$g * H * hx

  # Boundary condition: velocity at the calving front
  gamma <- (0.25 * p$A^(1/p$n) * (1 - p$rho/p$rhow) * p$rho * p$g * H[J+1])^p$n

  # Initial velocity (modify later to use ssainit())
  # u0 <- 0.001 * x
  # u0[1] <- u_bc # boundary condition at the ice divide (x1 = 0)
  # # u0 <- ssainit(p, x, beta, gamma, initchoice)
  u0 <- u

  Hstag <- 0.5 * (H[1:J] + H[2:(J+1)])
  tol <- 1.0e-14 # 1 / p$secpera
  eps_reg <- (1.0 / p$secpera) / p$L
  maxdiff <- Inf
  W <- numeric(J + 1)
  iter <- 0

  while (maxdiff > tol) {
    uxstag <- (u[2:(J+1)] - u[1:J]) / dx
    sqr_ux_reg <- uxstag^2 + eps_reg^2
    W[1:J] <- 2 * p$A^(-1/p$n) * Hstag * sqr_ux_reg^(((1/p$n) - 1) / 2)
    W[J+1] <- W[J]

    unew <- flowline(p$L, J, gamma, W, alpha, beta, u_bc, GL)
    maxdiff <- max(abs(unew - u))
    u <- unew
    iter <- iter + 1
    cat(".")

  }

  cat(sprintf("\nSSA solver did %d Picard iterations on dx = %.3f km grid\n", iter, dx / 1000))
  return(list(u = u, u0 = u0, zs = h, b = b, H = H, GL = GL))
}
