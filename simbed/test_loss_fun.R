## Test loss function

setwd("~/SSA_model/CNN/simbed/")

library(Matrix)
library(mvtnorm)
library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)


source("./source/custom_loss_function.R")

## Test loss function 

n_basis_funs <- 4L #dim(train_data$fric_coefs)[2]
n_gl <- 3L #im(train_data$grounding_line)[2]
d <- n_basis_funs * 2 + n_gl ## dimension of mean

Lb <- 2 * diag(n_basis_funs)
Lc <- diag(n_basis_funs)
Lg <- diag(n_gl)
# Lb[lower.tri(Lb, diag = F)] <- rep(0.5, n_basis_funs*(n_basis_funs-1)/2)
# Lc[lower.tri(Lc, diag = F)] <- rep(0.5, n_basis_funs*(n_basis_funs-1)/2)
# Lg[lower.tri(Lg, diag = F)] <- rep(0.5, n_gl*(n_gl-1)/2)

diag(Lb[-1, ]) <- rep(0.5, n_basis_funs-1)
diag(Lc[-1, ]) <- rep(0.5, n_basis_funs-1)
diag(Lg[-1, ]) <- rep(0.5, n_gl-1)

L <- bdiag(Lb, Lc, Lg)

y_true <- rep(1, d)
y_pred <- c(rep(0, d), log(diag(Lb)), diag(Lb[-1, ]), 
                        log(diag(Lc)), diag(Lc[-1, ]), 
                        log(diag(Lg)), diag(Lg[-1, ]))

test_loss <- posterior_loss_wrap(n_basis_funs, n_gl)(y_true, y_pred) 

Linv <- as.matrix(solve(L))
true_loss <- dmvnorm(y_true, mean = y_pred[1:d],
                    sigma = t(Linv) %*% Linv, log = T) + d/2 * log(2*pi)

## Tensorflow version

### Test construction of L matrix
# Lb_elems <- tf$constant(c(log(diag(Lb)), diag(Lb[-1, ])))
# Lb_tf <- construct_L_matrix_tf(Lb_elems, as.integer(n_basis_funs))
# Lb_tf

y_true_tf <- tf$constant(y_true)
y_pred_tf <- tf$constant(y_pred)

test_loss_tf <- posterior_loss_wrap_tf_og(n_basis_funs, n_gl)(y_true_tf, y_pred_tf) # - d/2 * log(2*pi)


# batch_size <- 2L
# y_true_tf2 <- tf$reshape(y_true_tf, shape = c(1L, dim(y_true_tf)))
# y_true_batch <- tf$tile(y_true_tf2, multiples = c(batch_size, 1L))
# y_pred_tf2 <- tf$reshape(y_pred_tf, shape = c(1L, dim(y_pred_tf)))
# y_pred_batch <- tf$tile(y_pred_tf2, multiples = c(batch_size, 1L))

# browser()

# test_loss_tf <- posterior_loss_wrap_tf(n_basis_funs, n_gl)(y_true_batch, y_pred_batch) # - d/2 * log(2*pi)

