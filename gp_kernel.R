# for matrix stuff

library('MASS')
library('emdbook')

dyn.load("kernel.so")
nr = 1e2
nc = 5
vec = rep(0, nr)
samples = mvrnorm(nr, c(rep(0, nc)), diag(nc))
mat = samples

f <- function(r, c, v, m){
    .C("kernel", row=as.integer(r), column=as.integer(c), mat=m, vec=v)
    }

out = f(nr, nc, vec, mat)

dmvnorm(samples, c(rep(0, nc)), diag(nc))

