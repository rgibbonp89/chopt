# for matrix stuff

library('MASS')
library('pracma')

dyn.load("kernel_rbf.so")
nrx = 10
nry = 20 
nc = 5
set.seed(10)
mat_x = mvrnorm(nrx, c(rep(0, nc)), diag(nc))
set.seed(10)
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
    rv_train = rep(0, xr*xr)
    rv_traintest = rep(0, xr*yr)
    rv_test = rep(0, yr*yr)
 
    out = .C("run", param=param, xrow=as.integer(xr), yrow=as.integer(yr), 
             column=as.integer(c), mat_x=mx, mat_y=my, 
             matmul_train=rv_train, 
             prod_train=rv_train, 
             kernel_train=rv_train,
             matmul_traintest=rv_traintest,
	     prod_traintest=rv_traintest,
	     kernel_traintest=rv_traintest,
             matmul_test=rv_test,
             prod_test=rv_test,
             kernel_test=rv_test)

    # reassign for output    
    out$matmul_train = matrix(out$matmul_train, nrow = xr, ncol = xr, byrow = T)
    out$prod_train = matrix(out$prod_train, nrow = xr, ncol = xr, byrow = T)
    out$kernel_train = matrix(out$kernel_train, nrow = xr, ncol = xr, byrow = T)
    
    out$matmul_traintest = matrix(out$matmul_traintest, nrow = xr, ncol = yr, byrow = T)
    out$prod_traintest = matrix(out$prod_traintest, nrow = xr, ncol = yr, byrow = T)
    out$kernel_traintest = matrix(out$kernel_traintest, nrow = xr, ncol = yr, byrow = T)

    out$matmul_test = matrix(out$matmul_test, nrow = yr, ncol = yr, byrow = T)
    out$prod_test = matrix(out$prod_test, nrow = yr, ncol = yr, byrow = T)
    out$kernel_test = matrix(out$kernel_test, nrow = yr, ncol = yr, byrow = T)

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
