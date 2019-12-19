# for matrix stuff

library('MASS')
library('emdbook')

dyn.load("kernel.so")
nr = 20
nc = 5
vec = rep(0, nr)
samples = mvrnorm(nr, c(rep(0, nc)), diag(nc))
mat = samples

mat = c(0,0, 1, 2)

f <- function(r, c, v, m){
    if(length(m) > 1){
    m <- c(t(m))
    }
    .C("kernel", row=as.integer(r), column=as.integer(c), mat=m, vec=v)
    }

out = f(nr, nc, vec, mat)






