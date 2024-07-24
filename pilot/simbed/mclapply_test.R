## Test mclapply error

library(parallel)

r <- mclapply(1:5, function(i) {
        if (i == 3L)
                stop("error in this process!")
        else
                return("success!")
}, mc.cores = 5)

str(r)

class(r[[3]])

inherits(r[[3]], "try-error")

bad <- sapply(r, inherits, what = "try-error")
bad

good <- r[!bad]
str(good)
