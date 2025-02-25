# posterior_loss_wrap_tf <- function {
  posterior_loss_tf <- function(y_true, y_pred) {
    ## Implement a loss of the form Gau(y_true; mean, cov)
    ## where mean and cov are constructed from y_pred
    # K <- backend() # call tensorflow
    # b <- y_true$shape[1] #dim(y_true)[1]
    # d <- as.integer(n_basis_funs * 2 + n_gl) ## Make this an argument in a wrapper function in the future
   
    # mean <- tf$constant(y_pred[, 1:d])
    # chol_vals <- tf$constant(y_pred[, (d+1):length(y_pred)])
    mean <- y_pred[, 1] ## this has dimension [batch_size, dim(output)]
    # chol_vals <- y_pred[, (d+1):dim(y_pred)[2]]
    sd <- exp(y_pred[, 2])
    
    mean <- tf$expand_dims(mean, axis = -1L)
    sd <- tf$expand_dims(sd, axis = -1L)

    # loss <- tf$math$square(y_true - mean)

    loss <- tf$math$square(y_true - mean) / (2 * tf$math$square(sd)) + 
             tf$math$log(sd) + 0.5 * tf$math$log(2 * tf$constant(pi))

    # loss <- tf$reduce_mean(loss, axis = 0L)

    # rhs <- tf$expand_dims(y_true - mean, axis = -1L)
    # v <- tf$matmul(tf$linalg$matrix_transpose(L_tf), rhs)
    # p1 <- tf$multiply(0.5, tf$matmul(tf$linalg$matrix_transpose(v), v))
    # p2 <- - tf$math$reduce_sum(tf$math$log(tf$linalg$diag_part(L_tf)), axis = 1L)
    
    # loss <- tf$squeeze(p1) + p2 #+ tf$math$multiply(0.5, tf$math$log(2 * tf$constant(pi)))
    return(loss)
  }
  # return(posterior_loss_tf)
# }


# true <- runif(3, -0.99, 0.99)#c(0.5, 0.3, 0.9)
# pred <- matrix(c(0.4, 0.2, 0.5, 0.1, 0.2, 0.3), nrow = 3)
# test <- lapply(1:3, function(i) dnorm(true[i], pred[i, 1], pred[i, 2], log = T))
# print(test)

# true_tf <- tf$constant(true, dtype = tf$float32)
# true_tf <- tf$expand_dims(true_tf, axis = -1L)
# pred_tf <- tf$constant(pred, dtype = tf$float32)
# print(posterior_loss_tf(true_tf, pred_tf))

