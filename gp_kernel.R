# for matrix stuff

library('MASS')
library('pracma')

dyn.load("kernel.so")
nrx = 20
nry = 10
nc = 5
mat_x = mvrnorm(nrx, c(rep(0, nc)), diag(nc))
mat_y = mvrnorm(nry, c(rep(0, nc)), diag(nc))


f <- function(mx, my){
    xr = nrow(mx)
    yr = nrow(my)
    c = ncol(mx)
    if(length(mx) > 1){
    mx <- c(t(mx))
    }
    if(length(my) > 1){
    my <- c(t(my))
    }
    rv = rep(0, xr*yr) 
    out = .C("rbf_kernel", xrow=as.integer(xr), yrow=as.integer(yr), 
    column=as.integer(c), mat_x=mx, mat_y=my, mat_res=rv)
    output = out$mat_res
    output_mat = matrix(output, nrow = xr, ncol = yr, byrow = T)
    return(output_mat)
    }

out = f(mat_x, mat_y)




