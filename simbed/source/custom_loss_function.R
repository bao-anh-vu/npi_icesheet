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

# quantile_loss_wrap <- function(tau) {
#   quantile_loss <- function(y_true, y_pred) {
#     K <- backend() ## this calls tensorflow

#     K$maximum(tau * (y_true - y_pred), (tau - 1) * (y_true - y_pred))
#   #   K$mean( K$abs( K$log( K$relu(y_true *1000 ) + 1 ) - 
#   #       K$log( K$relu(y_pred*1000 ) + 1)))
#   } 
#   return(quantile_loss) 
# }

quantile_loss_wrap <- function(tau) {
  quantile_loss <- function(y_true, y_pred) {
    K <- backend() ## this calls tensorflow
    bool_val <- K$greater(y_pred, y_true)
    ind_val <- tf$where(bool_val, 1, 0)

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
