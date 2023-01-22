# for matrix stuff

library('MASS')
library('pracma')

dyn.load("kernel.so")
nrx = 20
nry = 10 
nc = 5
mat_x = mvrnorm(nrx, c(rep(0, nc)), diag(nc))
mat_y = mvrnorm(nry, c(rep(0, nc)), diag(nc))


f <- function(mx, my, param){
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
    out = .C("rbf_kernel", param=param, xrow=as.integer(xr), yrow=as.integer(yr), 
    column=as.integer(c), mat_x=mx, mat_y=my, mat_res=rv, mat_prod=rv, mat_kernel=rv)
    out$mat_prod = matrix(out$mat_prod, nrow = xr, ncol = yr, byrow = T)
    out$mat_res = matrix(out$mat_res, nrow = xr, ncol = yr, byrow = T)
    out$mat_kernel = matrix(out$mat_kernel, nrow = xr, ncol = yr, byrow = T)
    return(out)
    }

start_time <- Sys.time()
out = f(mat_x, mat_y, 1.0)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
test <- mat_x%*%t(mat_y)
end_time <- Sys.time()
end_time - start_time



y = matrix(c(-0.23, -.21, 1.2, -0.23, 1.87, 0.3), nrow =2,  byrow=T)
x = matrix(c(0.01, 0.05, 0.12, 0.32, 1.20, -2.1), nrow = 2, byrow=T)
