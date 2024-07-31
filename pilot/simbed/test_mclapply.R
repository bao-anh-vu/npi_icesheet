library(parallel)

r <- mclapply(1:5, function(i) {
         if(i == 3L)
                stop("error in this process!")
        else
                return("success!")
}, mc.cores = 5)

## Check return value
str(r)

bad <- sapply(r, inherits, what = "try-error")
bad

r.good <- r[!bad]
str(r.good)

