# Model for predicting the mean and covariance of the posterior
create_model_posterior <- function(input_dim, output_dim, d) {
  model <- keras_model_sequential() %>%
    # layer_conv_1d(
    #     filters = 32, kernel_size = 20,
    #     padding = "same", activation = "relu",
    #     input_shape = input_dim #c(1000, 1)
    # ) %>%
    # layer_max_pooling_1d(pool_size = 2) %>%
    layer_conv_1d(
        filters = 64, kernel_size = 50,
        padding = "same", activation = "relu",
        input_shape = input_dim #c(1000, 1)
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_conv_1d(
        filters = 128, kernel_size = 20,
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_conv_1d(
        filters = 256, kernel_size = 10,
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_flatten() %>%
    # layer_dropout(rate = 0.1) %>%
    # layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = output_dim)

    model %>% compile(
        loss = posterior_loss_wrap_tf(d), #losses[1],
        optimizer = optimizer_adam(learning_rate = 0.00001),
        metrics = posterior_loss_wrap_tf(d) #c("mae") 
    )
    model
}

## Predicting just the mean
create_model <- function(input_dim, output_dim) {
  model <- keras_model_sequential() %>%
    layer_conv_1d(
        filters = 32, kernel_size = 5,
        padding = "same", activation = "relu",
        input_shape = input_dim #c(2001, 2, 2)
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_conv_1d(
        filters = 64, kernel_size = 5,
        padding = "same", activation = "relu",
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    layer_conv_1d(
        filters = 128, kernel_size = 3,
        padding = "same", activation = "relu"
    ) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
    # layer_conv_1d(
    #     filters = 256, kernel_size = 3,
    #     padding = "same", activation = "relu") %>%
    # layer_max_pooling_1d(pool_size = 2) %>%
    layer_flatten() %>%
    # layer_dropout(rate = 0.4) %>%
    # layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = output_dim)

    model %>% compile(
        loss = "mse",
        optimizer = optimizer_adam(learning_rate = 0.0001),
        metrics = list("mean_squared_error")
    )
    model
}