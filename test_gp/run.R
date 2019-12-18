dyn.load("./exe/gaussian_process.so")
x = rnorm(10, -3, 1)
l = length(x)
out_o <- .C("run", n=as.integer(l), x=x)


# for matrix stuff

dyn.load("gaussian_process.so")
nr = 5e3
nc = 5e4
vec = rep(0, nr)
mat = matrix(c(rnorm(nr*nc, -2, 1)), nrow = nr)

f <- function(r, c, v, m){
    .C("run", row=as.integer(r), column=as.integer(c), vec=v, mat=m)
    }

out = f(nr, nc, vec, mat)
