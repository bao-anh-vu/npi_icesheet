library(Matrix)  # for sparse matrices

flowline <- function(L, J, gamma, W, alpha, beta, u1, GL) {

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
  b[1] <- u1 #ug # boundary condition at the ice divide (x1 = 0)
  b[J+1] <- b[J+1] - 2 * gamma * dx * W[J+1]
  b[GL] <- 0 ## grounding line condition

  # # Initialize sparse matrix A (Bueler version)
  # A <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(J+1, J+1))

  # A[1,1] <- 1.0
  # for (j in 2:J) {
  #   A[j, j-1] <- W[j-1]
  #   A[j, j]   <- - (W[j-1] + W[j] + alpha[j] * dx^2)
  #   A[j, j+1] <- W[j]
  # }

  # A[J+1, J]   <- W[J] + W[J+1]
  # A[J+1, J+1] <- - (W[J] + W[J+1] + alpha[J+1] * dx^2)

  # # To make the velocity continuous across the grounding line,
  # # take u(xg + 1) as the average of u(xg + 2) and u(xg)
  # if (GL <= (J - 1)) {
  #   A[GL, (GL-1):(GL+1)] <- c(-1, 2, -1)
  # }

  # # Scale A by rows to improve conditioning
  # scale <- apply(abs(A), 1, max)
  # for (j in 1:(J+1)) {
  #   A[j, ] <- A[j, ] / scale[j]
  # }

#------------------
## Andrew's code for constructing the matrix A:

    v1 <- c(0, W[1 : J-1], W[J] + W[J + 1] ) # first element is dummy
    v2 <- c(0, -(W[1:(J-1)] + W[2:J] + alpha[2:J] * dx^2), -(W[J] + W[J + 1] + dx^2 * alpha[J + 1]))
    v3 <- c(0, W[2 : J], 0) # last element is dummy
    
    #A_mat[GL, (GL-1):(GL+1)] <- c(-1, 2, -1)
    v1[GL] <- -1
    v2[GL] <- 2
    v3[GL] <- -1

    #A_mat[1, 1] <- 1
    v2[1] <- 1
    
    V <- cbind(v1, v2, v3)
    
    # Row scaling to avoid singularity
    row_max <- as.vector(rowMaxs(abs(V))) #find maximum of each row
    V[1:nrow(V),] <- V[1:nrow(V),] / row_max[1:nrow(V)]
    
    A <- bandSparse(J + 1, J + 1, #dimensions
                        (-1):1, #band, diagonal is number 0
                        list(V[2:(J+1),1], V[,2], V[(1:J),3]))
#------------------

  # Scale b as well
  b <- b / row_max

  # Solve the linear system
  u <- solve(A, b)

  return(u)
}


