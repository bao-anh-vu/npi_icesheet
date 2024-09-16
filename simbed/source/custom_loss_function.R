# Custom loss function

# Mean Log Absolute Error
MLAE <- function( y_true, y_pred ) {
  K <- backend()
  K$mean( K$abs( K$log( K$relu(y_true *1000 ) + 1 ) - 
      K$log( K$relu(y_pred*1000 ) + 1)))
}

# Mean Squared Log Absolute Error
MSLAE <- function( y_true, y_pred ) {
  K <- backend()
  K$mean( K$pow( K$abs( K$log( K$relu(y_true *1000 ) + 1 ) - 
    K$log( K$relu(y_pred*1000 ) + 1)), 2))
}

## Quantile loss
quantile_loss <- function(y_true, y_pred, tau) {
  K <- backend() ## this calls tensorflow

    K$maximum(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
#   K$mean( K$abs( K$log( K$relu(y_true *1000 ) + 1 ) - 
#       K$log( K$relu(y_pred*1000 ) + 1)))
}


quantile_loss_wrap <- function(tau) {
  quantile_loss <- function(y_true, y_pred) {
    K <- backend() ## this calls tensorflow
    bool_val <- K$greater(y_pred, y_true)
    ind_val <- K$where(bool_val, 1, 0)

    # taus <- K$constant(rep(tau, dim(ind_val)), dtype = ind_val$dtype)
    loss <- (y_pred - y_true) * (ind_val - tau)
    # K$maximum(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
    return(loss)
  } 
  return(quantile_loss) 
}

## Compare pinball and quantile loss -- should be the same
# quantile_loss1 <- function(y_true, y_pred, tau) {

#     # K$maximum(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
#     max(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
    
# }

# quantile_loss2 <- function(y_true, y_pred, tau) {

#     # K$maximum(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
#     # max(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
#     (y_pred - y_true) * (as.numeric(y_true < y_pred) - tau)
# }

# ytrue <- 1 #rep(1, 10)
# ypred <- rnorm(1) #rnorm(10)

# q <- 0.95
# test1 <- quantile_loss1(ytrue, ypred, q)
# test2 <- quantile_loss2(ytrue, ypred, q)

#################
# Custom loss function for mean and covariance
posterior_loss_wrap <- function(n_basis_funs, n_gl) {
  posterior_loss <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred

    
    d <- n_basis_funs * 2 + n_gl ## Make this an argument in a wrapper function in the future
    
    mean <- y_pred[1:d]
    
    # n_cov_elements <- d + (150-1)*2 + (20-1)
    chol_vals <- y_pred[(d+1):length(y_pred)]

    ## Construct covariance matrix
    Lb_elems <- n_basis_funs + (n_basis_funs - 1)
    Lc_elems <- n_basis_funs + (n_basis_funs - 1)
    Lg_elems <- n_gl + (n_gl - 1)


    Lb <- construct_L_matrix(chol_vals[1:Lb_elems], n_basis_funs)
    Lc <- construct_L_matrix(chol_vals[(Lb_elems+1):(Lb_elems+Lc_elems)], n_basis_funs)
    Lg <- construct_L_matrix(chol_vals[(Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)

    L <- bdiag(Lb, Lc, Lg)

    v <- solve(solve(t(L)), y_true - mean)

    loss <- - 0.5 * t(v) %*% v + sum(log(diag(L)))
  }
  return(posterior_loss)
}


## Construct L matrix
construct_L_matrix <- function(vals, n) {
  L <- diag(exp(vals[1:n]))
  # L <- matrix(0, n, n)
  diag(L[-1, ]) <- vals[-(1:n)]
  # diag(L) <- exp(diag(L))
  return(L)
}

## Now do the same thing but in tensorflow

## Original version (not taking into account the batch size)
posterior_loss_wrap_tf_og <- function(n_basis_funs, n_gl) {
  posterior_loss_tf_og <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    # K <- backend() # call tensorflow
    d <- n_basis_funs * 2 + n_gl ## Make this an argument in a wrapper function in the future

    # mean <- tf$constant(y_pred[, 1:d])
    # chol_vals <- tf$constant(y_pred[, (d+1):length(y_pred)])
    mean <- y_pred[1:d] ## this has dimension [batch_size, dim(output)]
    chol_vals <- y_pred[(d+1):length(y_pred)]
    
    ## Construct covariance matrix
    Lb_elems <- n_basis_funs + (n_basis_funs - 1)
    Lc_elems <- n_basis_funs + (n_basis_funs - 1)
    Lg_elems <- n_gl + (n_gl - 1)

    Lb_tf <- construct_L_matrix_tf_og(chol_vals[1:Lb_elems], n_basis_funs)
    Lc_tf <- construct_L_matrix_tf_og(chol_vals[(Lb_elems+1):(Lb_elems+Lc_elems)], n_basis_funs)
    Lg_tf <- construct_L_matrix_tf_og(chol_vals[(Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)

    L_tf <- tf$linalg$LinearOperatorBlockDiag(list(Lb_tf, Lc_tf, Lg_tf))$to_dense()
    
    Linvt <- tf$linalg$inv(tf$transpose(L_tf))
    
    rhs <- tf$reshape(y_true - mean, c(as.integer(d), 1L))
    v <- tf$linalg$solve(Linvt, rhs)

    # loss <- - 0.5 * t(v) %*% v + sum(log(diag(L)))
    loss <- tf$multiply(-0.5, tf$matmul(tf$transpose(v), v)) + tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)))
  }
  return(posterior_loss_tf_og)
}

construct_L_matrix_tf_og <- function(vals, n) {
  L <- tf$linalg$diag(tf$exp(vals[1:n])) # diagonal matrix
  L2 <- tf$linalg$diag(vals[(n+1):length(vals)]) # smaller diag matrix, containing sub-diagonals
  MASK <- tf$constant(matrix(c(1L, 0L, 0L, 1L), 2L, 2L)) # padding for the smaller diag matrix so that it has the same dim as the bigger one

  L <- L + tf$pad(L2, MASK)
  L <- tf$linalg$LinearOperatorFullMatrix(L)
  return(L)
}

## Version that takes into account the batch size

posterior_loss_wrap_tf <- function(n_basis_funs, n_gl) {
  posterior_loss_tf <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    # K <- backend() # call tensorflow
    d <- n_basis_funs * 2 + n_gl ## Make this an argument in a wrapper function in the future
   
    # mean <- tf$constant(y_pred[, 1:d])
    # chol_vals <- tf$constant(y_pred[, (d+1):length(y_pred)])
    mean <- y_pred[, 1:d] ## this has dimension [batch_size, dim(output)]
    chol_vals <- y_pred[, (d+1):dim(y_pred)[2]]
    
    ## Construct covariance matrix
    Lb_elems <- n_basis_funs + (n_basis_funs - 1)
    Lc_elems <- n_basis_funs + (n_basis_funs - 1)
    Lg_elems <- n_gl + (n_gl - 1)

    Lb_tf <- construct_L_matrix_tf(chol_vals[, 1:Lb_elems], n_basis_funs)
    Lc_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+1):(Lb_elems+Lc_elems)], n_basis_funs)
    Lg_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)

    L_tf <- tf$linalg$LinearOperatorBlockDiag(list(Lb_tf, Lc_tf, Lg_tf))$to_dense()
    
    Linvt <- tf$linalg$inv(tf$transpose(L_tf))
    
    rhs <- tf$reshape(y_true - mean, c(as.integer(d), 1L))
    v <- tf$linalg$solve(Linvt, rhs)

    # loss <- - 0.5 * t(v) %*% v + sum(log(diag(L)))
    loss <- tf$multiply(-0.5, tf$matmul(tf$transpose(v), v)) + tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)))
  }
  return(posterior_loss_tf)
}

construct_L_matrix_tf <- function(vals, n) {
  b <- dim(vals)[1]
  L <- tf$linalg$diag(tf$exp(vals[, 1:n])) # diagonal matrix
  
  L2 <- tf$linalg$diag(vals[, (n+1):dim(vals)[2]]) # smaller diag matrix, containing sub-diagonals
  # MASK <- tf$constant(matrix(c(1L, 0L, 0L, 1L), 2L, 2L)) # padding for the smaller diag matrix so that it has the same dim as the bigger one
  MASK <- tf$constant(matrix(c(0L, 1L, 0L, 0L, 0L, 1L), 3L, 2L)) 
  # The 0-1-0-0-0-1 matrix is designed to pad each of 
  # the smaller diagonal matrices in the batch (L2) into 
  # the size of each of the the bigger diagonal matrices in L
  # but why it's 3x2 or why it's 0-1-0-0-0-1 is not clear to me
  # I just experimented with different values until it worked

  L <- L + tf$pad(L2, MASK)
  L <- tf$linalg$LinearOperatorFullMatrix(L)
  return(L)
}