solve_stress_balance <- function(L, J, gamma, W, alpha, beta, ug) {
  # % FLOWLINE  finite difference solution to 1D elliptic PDE problem (i.e.
  # % 2-point BVP)
  # %     (W(x) u_x)_x - alpha(x) u = beta(x)     on  0 < x < L
  # %     u(0) = ug,   u_x(L) = gamma
  # % form:
  # %   u = flowline(L,J,gamma,W,alpha,beta,ug)
  # % where:
  # %   L     = length of domain
  # %   J     = number of subintervals
  # %   gamma = right-hand Neumann boundary value
  # %   W     = vector of length J+1 (staggered grid values W_{j+1/2}, j=1:J+1)
  # %           (yes, last value should be at L+dx/2, past the end)
  # %   alpha = vector of length J+1 (regular grid values alpha_j, j=1:J+1)
  # %   beta  = same size as alpha
  # %   ug    = left-hand Dirichlet boundary value
  # % returns:
  # %   u     = solution vector (regular grid values u_j, j=1:J+1)
  # % examples (called by): TESTFLOWLINE, CONVANALYSIS, SSAFLOWLINE
  
  dx <- L / J
  b <- dx^2 * beta
  b[1] <- ug
  b[J+1] <- b[J+1] - 2 * gamma * dx * W[J+1]

  # Initialize sparse matrix A
  A <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(J+1, J+1))

  A[1,1] <- 1.0
  for (j in 2:J) {
    A[j, j-1] <- W[j-1]
    A[j, j]   <- - (W[j-1] + W[j] + alpha[j] * dx^2)
    A[j, j+1] <- W[j]
  }

  A[J+1, J]   <- W[J] + W[J+1]
  A[J+1, J+1] <- - (W[J] + W[J+1] + alpha[J+1] * dx^2)

  # Scale A by rows to improve conditioning
  scale <- apply(abs(A), 1, max)
  for (j in 1:(J+1)) {
    A[j, ] <- A[j, ] / scale[j]
  }
  b <- b / scale

  # Solve the linear system
  u <- solve(A, b)
  return(u)
}
