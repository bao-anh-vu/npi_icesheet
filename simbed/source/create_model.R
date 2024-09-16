# Define a simple sequential model
create_model <- function(input_dim, output_dim) {
  model <- keras_model_sequential() %>%
    layer_conv_2d(
        filters = 32, kernel_size = c(5, 5),
        padding = "same", activation = "relu",
        input_shape = input_dim #c(2001, 2, 2)
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 64, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 128, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    # layer_conv_2d(filters = 256, kernel_size = c(3, 3),
    #                 padding = "same", activation = "relu") %>%
    # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_flatten() %>%
    layer_dropout(rate = 0.5) %>%
    # layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = output_dim)

    model %>% compile(
        loss = "mse",
        optimizer = optimizer_adam(learning_rate = 0.001),
        metrics = list("mean_squared_error")
    )
    model
}

create_model_quantile <- function(input_dim, output_dim, quantile) {
  model <- keras_model_sequential() %>%
    layer_conv_2d(
        filters = 32, kernel_size = c(5, 5),
        padding = "same", activation = "relu",
        input_shape = input_dim #c(2001, 2, 2)
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 64, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 128, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    # layer_conv_2d(filters = 256, kernel_size = c(3, 3),
    #                 padding = "same", activation = "relu") %>%
    # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_flatten() %>%
    layer_dropout(rate = 0.5) %>%
    # layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = output_dim)

    # Compile the model, and select one of the loss functions
    # q <- 0.95 #c(0.1, 0.5, 0.9)

    # losses <- c(keras::loss_mean_squared_error,  
    #     keras::loss_mean_squared_logarithmic_error, 
    #     MLAE, MSLAE, quantile_loss_wrap)

    # metrics <- c("mae", "mse", "msle", "mlae", "mslae")    
    model %>% compile(
        # optimizer = "rmsprop",
        loss = quantile_loss_wrap(quantile), #losses[1],
        optimizer = optimizer_adam(learning_rate = 0.0005),
        metrics = quantile_loss_wrap(quantile) #c("mae") # is this the right thing to do?
    )

    # model %>% compile(
    #     loss = "mse",
    #     optimizer = optimizer_adam(learning_rate = 0.001),
    #     metrics = list("mean_squared_error")
    # )
    model
}

create_model_posterior <- function(input_dim, output_dim, n_basis_funs, n_gl) {
  model <- keras_model_sequential() %>%
    layer_conv_2d(
        filters = 32, kernel_size = c(5, 5),
        padding = "same", activation = "relu",
        input_shape = input_dim #c(2001, 2, 2)
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 64, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_conv_2d(
        filters = 128, kernel_size = c(3, 3),
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    # layer_conv_2d(filters = 256, kernel_size = c(3, 3),
    #                 padding = "same", activation = "relu") %>%
    # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
    layer_flatten() %>%
    layer_dropout(rate = 0.5) %>%
    # layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = output_dim)

    model %>% compile(
        loss = posterior_loss_wrap_tf(n_basis_funs, n_gl), #losses[1],
        optimizer = optimizer_adam(learning_rate = 0.0005),
        metrics = posterior_loss_wrap_tf(n_basis_funs, n_gl) #c("mae") 
    )
    model
}

# # Define a simple sequential model
# create_model_bed <- function(input_dim, output_dim) {
#   model <- keras_model_sequential() %>%
#     layer_conv_2d(
#         filters = 32, kernel_size = c(5, 5),
#         padding = "same", activation = "relu",
#         input_shape = input_dim #c(2001, 50, 2)
#     ) %>%
#     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#     layer_conv_2d(
#         filters = 64, kernel_size = c(3, 3),
#         padding = "same", activation = "relu"
#     ) %>%
#     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#     layer_conv_2d(
#         filters = 128, kernel_size = c(3, 3),
#         padding = "same", activation = "relu"
#     ) %>%
#     layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#     # layer_conv_2d(filters = 256, kernel_size = c(3, 3),
#     #                 padding = "same", activation = "relu") %>%
#     # layer_max_pooling_2d(pool_size = c(2, 2)) %>%
#     layer_flatten() %>%
#     layer_dropout(rate = 0.25) %>%
#     # layer_dense(units = 256, activation = "relu") %>%
#     layer_dense(units = output_dim)

#     model %>% compile(
#         loss = "mse",
#         optimizer = optimizer_adam(learning_rate = 0.0005),
#         metrics = list("mean_squared_error")
#     )
#     model
# }

