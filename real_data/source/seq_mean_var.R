compute_mean_seq <- function(x) {
    if (!is.null(dim(x))) {
        x <- as.vector(x)
    }

    x_mean <- x[1]

    for (i in 2:length(x)) {
        x_mean <- x_mean + (x[i] - x_mean)/(i)
    }

    return(x_mean)
}

compute_var_seq <- function(x) {
    sum_of_sq <- 0
    x_var <- 0
    x_mean_prev <- x[1]
    
    for (i in 2:length(x)) {
        x_mean_curr <- x_mean_prev + (x[i] - x_mean_prev)/(i)
        sum_of_sq <- sum_of_sq + (x[i] - x_mean_prev)*(x[i] - x_mean_curr)
        x_mean_prev <- x_mean_curr
    }
    x_var <- sum_of_sq/(length(x)-1)
    return(x_var)
}