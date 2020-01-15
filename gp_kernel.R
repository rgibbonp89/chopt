library('MASS')
library('pracma')
library('e1071')
library('caret')
library('mixtools')


# Load compiled C
dyn.load("kernel_rbf.so")

# Computes mu and vcov matrix for given test set
#' @param mx training set (without y)
#' @param my test set (without y)
#' @param y test set y
f <- function(mx, my, y, param){

    # Calculate dimensions
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

    # Flatten matrices to vectors
    rv_train = rep(0, xr*xr)
    rv_traintest = rep(0, xr*yr)
    rv_test = rep(0, yr*yr)
    mu_out = rep(0, yr)
    vcov_out = rep(0, yr*yr)
 
    # Run main C function - annotation can be found in .c file itself
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

    return(list(mu = out$mu_out, vcov = out$vcov_out))
    }



#' Fit SVM to training data (only iris for the time being) and predict from validation set
#' @param train_data training set (including target)
#' @param validation_data validation set (iincluding target)
#' @param C C parameter of SVM
#' @param gamma gamma parameter of SVM
svm_fit_predict <- function(train_data, validation_data, C, gamma){

   # At this stage, just use simple performance metric to validate

    model <- svm(Species ~ ., data=train_data, gamma=gamma, cost=C)
    validation_x <- subset(validation_data, select=-Species)
    validation_y <- subset(validation_data, select=Species)
    predictions <- predict(model, validation_x)

    return(sum(validation_y$Species == predictions)/length(predictions))
}



#' Optimize SVM hyperparameters (IN PROGRESS)
#' Gamma parameter - polynomial order (essentially for linear decision boundary)
#' C parameter - penalty for misclassifying

#' @param train training set (including target)
#' @param validation validation set (iincluding target)
#' @param C_lower lower bound for C hyperparameter of SVM
#' @param C_upper upper bound for C hyperparameter of SVM
#' @param C_num number of candidates to try for C hyperparameter of SVM
#' @param gamma_lower lower bound for gamma hyperparameter of SVM
#' @param gamma_upper upper bound for gamma hyperparameter of SVM
#' @param gamma_num number of candidates to try for gamma hyperparameter of SVM
#' @param warmups number of warmup steps to run
#' @param iterations total number of iterations to run
optimize_params_svm <- function(train, validation, C_lower, C_upper, C_num, gamma_lower, gamma_upper, gamma_num, warmups=5, iterations = 5){


   C_candidates <- sort(runif(C_num, C_lower, C_upper))
   gamma_candidates <- sort(runif(gamma_num, gamma_lower, gamma_upper))
   C_gamma_combi <- expand.grid(list(C_candidates=C_candidates, gamma_candidates=gamma_candidates))

   # Run warmups (fixed at 5 for now)

   inputs <- c()
   score <- c()
   for(warmup in 1:warmups){
   
       index <- sample(dim(C_gamma_combi)[1], 1)
       cg_tmp <- C_gamma_combi[index,]
       acc_tmp <- svm_fit_predict(train, validation, cg_tmp$C_candidates, cg_tmp$gamma_candidates)
       inputs <- rbind(inputs, cg_tmp)
       score <- rbind(score, acc_tmp)
       C_gamma_combi <- C_gamma_combi[-index,]
 
   }
   

   out <- f(inputs, C_gamma_combi, score, 10.0)
   best_ <- score[which.max(score),]
   next_move_cands <- expected_improvement(out, best_)
   
   for(iter in 1:iterations){
   
   cg_tmp <- C_gamma_combi[which.max(next_move_cands),]
   acc_tmp <- svm_fit_predict(train, validation, cg_tmp$C_candidates, cg_tmp$gamma_candidates)
   inputs <- rbind(inputs, cg_tmp)
   score <- rbind(score, acc_tmp)
   best_ <- score[which.max(score),]   
   out <- f(inputs, C_gamma_combi, score, 10.0)
   next_move_cands <- expected_improvement(out, best_)
   }

   return(list(scores=score, inputs=inputs))

}

#' Expected improvement for set of candidate hyperparameters
#' @param gp_res output from C function Gaussian Process
#' @param best best hyperparameter score so far
expected_improvement <- function(gp_res, best){
    diff <- best - gp_res$mu
    sigma <- diag(gp_res$vcov)
    return(diff*dnorm(diff/sigma, 0, 1) + sigma*pnorm(diff/sigma, 0, 1))
}



###########################
# Running a simple example
###########################


# Validation of above on IRIS data
all_data = subset(iris, Species == 'virginica' | Species == 'versicolor')

# TVT split
trainIndex <- createDataPartition(all_data$Species, p = .7,
                                  list = FALSE,
                                  times = 1)
train_set <- all_data[trainIndex,]
val_test <- all_data[-trainIndex,]

# Specify hyperparameter ranges
C_lower <- 1e-4
C_upper <- 10
C_num <- 20

gamma_lower <- 1e-4
gamma_upper <- 10
gamma_num <- 20

# Run the algorithm
next_move <- optimize_params_svm(train_set, val_test, C_lower, C_upper, C_num, gamma_lower, gamma_upper, gamma_num, iterations=30)



