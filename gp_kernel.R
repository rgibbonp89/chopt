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

betas = rnorm(nc)
y = mat_x%*%c(betas)


f <- function(mx, my, y, param){
    xr = nrow(mx)
    yr = nrow(my)
    c = ncol(mx)
    if(length(mx) > 1){
    mx <- c(t(mx))
    }
    if(length(my) > 1){
    my <- c(t(my))
    }
    if(length(y) > 1){
    y <- c(t(y))
    }

    rv_train = rep(0, xr*xr)
    rv_traintest = rep(0, xr*yr)
    rv_test = rep(0, yr*yr)
    mu_out = rep(0, yr)
    vcov_out = rep(0, yr*yr)
 
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
             kernel_test=rv_test,
             y_in=y,
             mu_out=mu_out,
             vcov_out=vcov_out)

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
    
    out$mu_out = matrix(out$mu_out, nrow = yr, ncol = 1, byrow = T)
    out$vcov_out = matrix(out$vcov_out, nrow = yr, ncol = yr, byrow = T)

    return(out)
    }

out = f(mat_x, mat_y, y, 1.0)

# mus match exactly!
K_inv = solve(out$kernel_train)
K_s = out$kernel_traintest
K_ss = out$kernel_test

mu = (t(K_s)%*%K_inv)%*%y
vcov = K_ss - (t(K_s)%*%K_inv)%*%K_s
