# Custom loss function for mean and covariance
posterior_loss_wrap <- function(n_bed_basis, n_fric_basis, n_gl) {
  posterior_loss <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred

    
    d <- n_bed_basis + n_fric_basis + n_gl ## Make this an argument in a wrapper function in the future
    
    mean <- y_pred[1:d]
    
    # n_cov_elements <- d + (150-1)*2 + (20-1)
    chol_vals <- y_pred[(d+1):length(y_pred)]

    ## Construct covariance matrix
    Lb_elems <- n_bed_basis + (n_bed_basis - 1)
    Lc_elems <- n_fric_basis + (n_fric_basis - 1)
    Lg_elems <- n_gl + (n_gl - 1)

    Lb <- construct_L_matrix(chol_vals[1:Lb_elems], n_bed_basis)
    Lc <- construct_L_matrix(chol_vals[(Lb_elems+1):(Lb_elems+Lc_elems)], n_fric_basis)
    Lg <- construct_L_matrix(chol_vals[(Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)

    L <- bdiag(Lb, Lc, Lg)

    v <- t(L) %*% (y_true - mean)

    loss <- 0.5 * t(v) %*% v - sum(log(diag(L)))
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
posterior_loss_wrap_tf_og <- function(n_bed_basis, n_fric_basis, n_gl) {
  posterior_loss_tf_og <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    # K <- backend() # call tensorflow
    d <- n_bed_basis + n_fric_basis + n_gl ## Make this an argument in a wrapper function in the future

    # mean <- tf$constant(y_pred[, 1:d])
    # chol_vals <- tf$constant(y_pred[, (d+1):length(y_pred)])
    mean <- y_pred[1:d] ## this has dimension [batch_size, dim(output)]
    chol_vals <- y_pred[(d+1):length(y_pred)]
    
    ## Construct covariance matrix
    Lb_elems <- n_bed_basis + (n_bed_basis - 1)
    Lc_elems <- n_fric_basis + (n_fric_basis - 1)
    Lg_elems <- n_gl + (n_gl - 1)

    Lb_tf <- construct_L_matrix_tf_og(chol_vals[1:Lb_elems], n_bed_basis)
    Lc_tf <- construct_L_matrix_tf_og(chol_vals[(Lb_elems+1):(Lb_elems+Lc_elems)], n_fric_basis)
    Lg_tf <- construct_L_matrix_tf_og(chol_vals[(Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)

    L_tf <- tf$linalg$LinearOperatorBlockDiag(list(Lb_tf, Lc_tf, Lg_tf))$to_dense()
    
    # Linvt <- tf$linalg$inv(tf$transpose(L_tf))
    
    rhs <- tf$reshape(y_true - mean, c(as.integer(d), 1L))
    # v <- tf$linalg$solve(Linvt, rhs)
    v <- tf$matmul(tf$linalg$matrix_transpose(L_tf), rhs)

    # loss <- - 0.5 * t(v) %*% v + sum(log(diag(L)))
    loss <- tf$multiply(0.5, tf$matmul(tf$transpose(v), v)) - tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)))
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

posterior_loss_wrap_tf <- function(n_bed_basis, n_fric_basis, n_gl) {
  posterior_loss_tf <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    # K <- backend() # call tensorflow
    # b <- y_true$shape[1] #dim(y_true)[1]
    d <- as.integer(n_bed_basis + n_fric_basis + n_gl) ## Make this an argument in a wrapper function in the future
   
    # mean <- tf$constant(y_pred[, 1:d])
    # chol_vals <- tf$constant(y_pred[, (d+1):length(y_pred)])
    mean <- y_pred[, 1:d] ## this has dimension [batch_size, dim(output)]
    chol_vals <- y_pred[, (d+1):dim(y_pred)[2]]
    
    ## Construct covariance matrix
    Lb_elems <- n_bed_basis + (n_bed_basis - 1)
    Lc_elems <- n_fric_basis + (n_fric_basis - 1)
    Lg_elems <- n_gl + (n_gl - 1)

    Lb_tf <- construct_L_matrix_tf(chol_vals[, 1:Lb_elems], n_bed_basis)
    Lc_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+1):(Lb_elems+Lc_elems)], n_fric_basis)
    Lg_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)
    
    L_tf <- tf$linalg$LinearOperatorBlockDiag(list(Lb_tf, Lc_tf, Lg_tf))$to_dense()
    
    # Linvt <- tf$linalg$inv(tf$linalg$matrix_transpose(L_tf))
    # rhs <- tf$reshape(y_true - mean, c(tf$shape(y_true), 1L)) # why is y_true of length 1280 not 20480 in this case?
    rhs <- tf$expand_dims(y_true - mean, axis = -1L)
    
    # v <- tf$linalg$solve(Linvt, rhs)

    v <- tf$matmul(tf$linalg$matrix_transpose(L_tf), rhs)

    p1 <- tf$multiply(0.5, tf$matmul(tf$linalg$matrix_transpose(v), v))
    # loss <- - 0.5 * t(v) %*% v + sum(log(diag(L)))
    p2 <- - tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)), axis = 1L)
    # loss <- p1 + tf$reshape(p2, c(dim(p2), 1L, 1L))
    loss <- tf$squeeze(p1) + p2
    
  }
  return(posterior_loss_tf)
}

construct_L_matrix_tf <- function(vals, n) {
  # b <- dim(vals)[1]
  L <- tf$linalg$diag(tf$exp(vals[, 1:n])) # diagonal matrix
  
  L2 <- tf$linalg$diag(vals[, (n+1):dim(vals)[2]]) # smaller diag matrix, containing sub-diagonals
  # MASK <- tf$constant(matrix(c(1L, 0L, 0L, 1L), 2L, 2L)) # padding for the smaller diag matrix so that it has the same dim as the bigger one
  MASK <- tf$constant(matrix(c(0L, 1L, 0L, 0L, 0L, 1L), 3L, 2L)) 
  # The 0-1-0-0-0-1 matrix is designed to pad each of 
  # the smaller diagonal matrices in the batch (L2) into 
  # the size of each of the bigger diagonal matrices in L
  # but why it's 3x2 or why it's 0-1-0-0-0-1 is not clear to me
  # I just experimented with different values until it worked

  L <- L + tf$pad(L2, MASK)
  L <- tf$linalg$LinearOperatorLowerTriangular(L)
  return(L)
}