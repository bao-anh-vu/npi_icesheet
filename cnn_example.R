## Example from ISLR
## Predict baseball players' salaries using a neural network

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ISLR)
library(ggplot2)

Gitters <- na.omit(Hitters)
n <- nrow(Gitters)
set.seed(13)
ntest <- trunc(n / 3)
testid <- sample(1:n, ntest)

x <- scale(model.matrix(Salary ~ . - 1, data = Gitters))
y <- Gitters$Salary

modnn <- keras_model_sequential() %>%
            layer_dense(units = 50, activation = "relu",
                        input_shape = ncol(x)) %>%
            layer_dropout(rate = 0.4) %>%
            layer_dense(units = 1)

modnn %>% compile(loss = "mse",
                    optimizer = optimizer_rmsprop(),
                    metrics = list("mean_absolute_error")
                    )

history <- modnn %>% fit(
                        x[-testid, ], y[-testid], epochs = 1500, batch_size = 32,
                        validation_data = list(x[testid, ], y[testid])
                        )
plot(history)

## Predict on the test data
npred <- predict(modnn, x[testid, ])
mean(abs(y[testid] - npred))
