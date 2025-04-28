ssaflowline <- function(p, J, H, b, ug, initchoice) {

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
  # %   ug = velocity boundary value at x=0 end of domain
  # %   initchoice = 1,2 chooses initial guess scheme;
  # %                use 1 for ice shelves and 2 for ice streams
  # % outputs:
  # %   u = the velocity solution (m s-1), a length J+1 column vector
  # %   u0 = the velocity initial guess (m s-1), same size as u
  # % Does Picard iteration to solve nonlinear SSA problem.
  # % Calls:  FLOWLINE to solve inner linear PDE boundary value problem,
  # %         SSAINIT  to get initial iterate

  # Check if all required arguments are provided
  if (missing(p) || missing(J) || missing(H) || missing(b) || missing(ug) || missing(initchoice)) {
    stop("All 6 input arguments are required.")
  }

  dx <- p$L / J
  x <- seq(0, p$L, length.out = J + 1)
  xstag <- seq(dx/2, p$L + dx/2, by = dx)

  alpha <- rep(p$C, length(x))  # Drag coefficient
  h <- H + b # surface elevation (z_s in my notation)

  # Compute regular grid values of slope
  hx <- numeric(length(H))
  hx[1] <- (h[2] - h[1]) / dx
  hx[2:J] <- (h[3:(J+1)] - h[1:(J-1)]) / (2 * dx)
  hx[J+1] <- (h[J+1] - h[J]) / dx

  beta <- p$rho * p$g * H * hx

  gamma <- (0.25 * p$A^(1/p$n) * (1 - p$rho/p$rhow) * p$rho * p$g * H[J+1])^p$n

  u0 <- ssainit(p, x, ug, beta, gamma, initchoice)
  u <- u0

  Hstag <- 0.5 * (H[1:J] + H[2:(J+1)])
  tol <- 1.0e-14 # 1 m/yr
  eps_reg <- (1.0 / p$secpera) / p$L
  maxdiff <- Inf
  W <- numeric(J + 1)
  iter <- 0

  u_list <- list()
  while (maxdiff > tol) {
    uxstag <- (u[2:(J+1)] - u[1:J]) / dx
    sqr_ux_reg <- uxstag^2 + eps_reg^2
    W[1:J] <- 2 * p$A^(-1/p$n) * Hstag * sqr_ux_reg^(((1/p$n) - 1) / 2)
    W[J+1] <- W[J]

    unew <- flowline(p$L, J, gamma, W, alpha, beta, ug)
    maxdiff <- max(abs(unew - u))
    u <- unew
    u_list[[iter + 1]] <- u
    iter <- iter + 1
    cat(".")
    browser()
  }

  cat(sprintf("\nSSA solver did %d Picard iterations on dx = %.3f km grid\n", iter, dx / 1000))
  return(list(u_list = u_list, u = u, u0 = u0, zs = h))
}
