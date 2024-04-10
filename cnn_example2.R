## Example from ISLR
## MNIST dataset

library(keras)
reticulate::use_condaenv("myenv", required = TRUE)
library(tensorflow)
library(ISLR)
library(ggplot2)

cifar100 <- dataset_cifar100()
names(cifar100)
x_train <- cifar100$train$x[1:1000,,,]
g_train <- cifar100$train$y[1:1000,]
x_test <- cifar100$test$x[1:100,,,]
g_test <- cifar100$test$y[1:100,]
dim(x_train)
range(x_train[1,,, 1])

x_train <- x_train / 255 ## standardize
x_test <- x_test / 255
y_train <- to_categorical(g_train, 100) ## one hot coding
dim(y_train)

# library(jpeg)
# par(mar = c(0, 0, 0, 0), mfrow = c(5, 5))
# index <- sample(seq(50000), 5)
# for (i in index) plot(as.raster(x_train[i,,, ]))

model <- keras_model_sequential() %>%
            layer_conv_2d(filters = 32, kernel_size = c(3, 3),
                            padding = "same", activation = "relu",
                            input_shape = c(32, 32, 3)) %>%
            layer_max_pooling_2d(pool_size = c(2, 2)) %>%
            layer_conv_2d(filters = 64, kernel_size = c(3, 3),
                            padding = "same", activation = "relu") %>%
            layer_max_pooling_2d(pool_size = c(2, 2)) %>%
            layer_conv_2d(filters = 128, kernel_size = c(3, 3),
                            padding = "same", activation = "relu") %>%
            layer_max_pooling_2d(pool_size = c(2, 2)) %>%
            layer_conv_2d(filters = 256, kernel_size = c(3, 3),
                            padding = "same", activation = "relu") %>%
            layer_max_pooling_2d(pool_size = c(2, 2)) %>%
            layer_flatten() %>%
            layer_dropout(rate = 0.5) %>%
            layer_dense(units = 512, activation = "relu") %>%
            layer_dense(units = 100, activation = "softmax")

summary(model)

model %>% compile(loss = "categorical_crossentropy",
                    optimizer = optimizer_rmsprop(), metrics = c("accuracy"))
history <- model %>% fit(x_train, y_train, epochs = 30,
                        batch_size = 128, validation_split = 0.2)
accuracy <- function(pred, truth)
    mean(drop(as.numeric(pred)) == drop(truth))

model %>% predict(x_test) %>% k_argmax() %>% accuracy(g_test)
