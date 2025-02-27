## TF version that takes into account the batch size

posterior_loss_wrap_tf <- function(d) {
  posterior_loss_tf <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    mean <- y_pred[, 1:d] ## this has dimension [batch_size, dim(output)]
    chol_vals <- y_pred[, (d+1):dim(y_pred)[2]]
    
    ## Construct covariance matrix
    # Lb_elems <- n_basis_funs + (n_basis_funs - 1)
    # Lc_elems <- n_basis_funs + (n_basis_funs - 1)
    # Lg_elems <- n_gl + (n_gl - 1)

    # Lb_tf <- construct_L_matrix_tf(chol_vals[, 1:Lb_elems], n_basis_funs)
    # Lc_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+1):(Lb_elems+Lc_elems)], n_basis_funs)
    # Lg_tf <- construct_L_matrix_tf(chol_vals[, (Lb_elems+Lc_elems+1):(Lb_elems+Lc_elems+Lg_elems)], n_gl)
    
    # L_tf <- tf$linalg$LinearOperatorBlockDiag(list(Lb_tf, Lc_tf, Lg_tf))$to_dense()
    
    L_tf <- construct_L_matrix_tf(chol_vals, d)
    
    rhs <- tf$expand_dims(y_true - mean, axis = -1L)
    v <- tf$matmul(tf$linalg$matrix_transpose(L_tf), rhs)
    p1 <- tf$multiply(0.5, tf$matmul(tf$linalg$matrix_transpose(v), v))
    p2 <- - tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)), axis = 1L)
    
    loss <- tf$squeeze(p1) + p2 #+ tf$math$multiply(0.5, tf$math$log(2 * tf$constant(pi)))
    loss <- tf$expand_dims(loss, axis = -1L)

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
  
  # L <- tf$linalg$LinearOperatorLowerTriangular(L)

  
  return(L)
}

## Non-TF version

# Custom loss function for mean and covariance
posterior_loss_wrap <- function(d) {
  posterior_loss <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    
    mean <- y_pred[1:d]
    chol_vals <- y_pred[(d+1):length(y_pred)]

    ## Construct covariance matrix
    L <- construct_L_matrix(chol_vals, d)

    ## Calculate loss function
    v <- t(L) %*% (y_true - mean)
    loss <- 0.5 * t(v) %*% v - sum(log(diag(L)))
  }
  return(posterior_loss)
}


## Construct L matrix
construct_L_matrix <- function(vals, n) {
  L <- diag(exp(vals[1:n]))
  # L <- matrix(0, n, n)
  L[lower.tri(L, diag = F)] <- vals[-(1:n)]
  # diag(L) <- exp(diag(L))
  return(L)
}

# ## Test loss function
# batchsize <- 3L
# d <- 2L
# y_true <- matrix(rnorm(batchsize * d), nrow = batchsize)
# y_pred <- matrix(rnorm(batchsize * (d + d + (d - 1))), nrow = batchsize) # d elements for mean, d + (d - 1) for chol_vals
# loss_val <- lapply(1:batchsize, function(i) posterior_loss_wrap(d)(y_true[i, ], y_pred[i, ]))
# print(loss_val)

# ## Loss from tf
# y_pred_tf <- tf$constant(y_pred, dtype = tf$float64)
# y_true_tf <- tf$constant(y_true, dtype = tf$float64)

# loss_tf <- posterior_loss_wrap_tf(d)(y_true = y_true_tf, y_pred = y_pred_tf)
# print(loss_tf)

