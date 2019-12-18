dyn.load("kernel.so")
nr = 5e2
nc = 5e1
vec = rep(0, nr)
mat = matrix(c(rnorm(nr*nc, -2, 1)), nrow = nr)

f <- function(r, c, v, m){
    .C("run", row=as.integer(r), column=as.integer(c), vec=v, mat=m)
    }

out = f(nr, nc, vec, mat)

# for matrix stuff

dyn.load("kernel.so")
nr = 1e2
nc = 5
vec = rep(0, nr)
mat = matrix(c(rnorm(nr*nc, 0, 1)), nrow = nr)

f <- function(r, c, v, m){
    .C("kernel", row=as.integer(r), column=as.integer(c), mat=m, vec=v)
    }

out = f(nr, nc, vec, mat)
