/*

RBF KERNEL FOR HYPERARAMETER TUNING

Designed to integrate with R through the .C function, this routine will calculate a mu vector and variance covariance matrix
for unseen hyperparameter combinations.

*/

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>

/*
@param_in = tuning parameter for RBF Kernel
@xrow = number of rows in xmatrix
@yrow = number of rows in ymatrix
@col = number of columns in both matrices (has to be the same in both, obviously)
@input_matrix_x = x matrix (training set)
@input_matrix_y = y matrix (test set)
@matmul_matrix_train, @product_matrix_train = (xrow, xrow) blank matrix to store intermediate results during matrix operations
@kernel_matrix_train = (xrow, xrow) matrix that will hold output of kernel calculation

Same goes for _traintest and _test

@y_vector = (xrow,) vector of train targets
@mu_vector = (yrow,) length vector of outputted mu estimates
@vcov_matrix_out = (yrow, yrow) estimated variance covariance matrix
*/

int run(double *param_in, int *xrow, int *yrow, int *col,
        double *input_matrix_x, double *input_matrix_y,
        double *matmul_matrix_train,
        double *product_matrix_train,
        double *kernel_matrix_train,
        double *matmul_matrix_traintest,
        double *product_matrix_traintest,
        double *kernel_matrix_traintest,
        double *matmul_matrix_test,
        double *product_matrix_test,
        double *kernel_matrix_test,
        double *y_vector,
        double *mu_vector_out,
        double *vcov_matrix_out
        )
{
    // Allocate matrices for storing output
    size_t xr = xrow[0];
    size_t yr = yrow[0];
    gsl_matrix * Ktrain = gsl_matrix_calloc(xr, xr);
    gsl_matrix * Ktraintest = gsl_matrix_calloc(xr, yr);
    gsl_matrix * Ktest = gsl_matrix_calloc(yr, yr);
    gsl_matrix * Kinv = gsl_matrix_calloc(xr, xr);

    // Pass pointers themselves as arguments to rbf_kernel()
    rbf_kernel(param_in, xrow, xrow, col, input_matrix_x,
               input_matrix_x, matmul_matrix_train, product_matrix_train,
               kernel_matrix_train);
    rbf_kernel(param_in, xrow, yrow, col, input_matrix_x,
               input_matrix_y, matmul_matrix_traintest, product_matrix_traintest,
               kernel_matrix_traintest);
    rbf_kernel(param_in, yrow, yrow, col, input_matrix_y,
               input_matrix_y, matmul_matrix_test, product_matrix_test,
               kernel_matrix_test);

    // Assign to matrices
    matrix_assign_GSL(xrow, xrow, kernel_matrix_train, Ktrain);
    matrix_assign_GSL(xrow, xrow, kernel_matrix_train, Kinv);
    matrix_assign_GSL(xrow, yrow, kernel_matrix_traintest, Ktraintest);
    matrix_assign_GSL(yrow, yrow, kernel_matrix_test, Ktest);

    // Calculate Ktrain^{-1}
    matrix_inverse(Kinv, Kinv, xr);

    // Calculate mu vector
    gsl_matrix * K_traintest_Kinv = gsl_matrix_calloc(yr, xr);
    gsl_matrix * y_vector_in = gsl_matrix_calloc(xr, 1);
    gsl_matrix * mu_vector = gsl_matrix_calloc(yr, 1);

    int i=0;
    for(i=0;i<xr;i++)
    {
        gsl_matrix_set(y_vector_in, i, 0, y_vector[i]);
    }

    gsl_blas_dgemm(CblasTrans, CblasNoTrans,
                  1.0, Ktraintest, Kinv,
                  0.0, K_traintest_Kinv);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, K_traintest_Kinv, y_vector_in,
                  0.0, mu_vector);

    int j=0;
    for(j=0;j<yr;j++)
    {
        mu_vector_out[j] = gsl_matrix_get(mu_vector, j, 0);
    }


    // Calculate covariance matrix
    gsl_matrix * vcovmatrix_holder = gsl_matrix_calloc(yr, yr);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                  1.0, K_traintest_Kinv, Ktraintest,
                  0.0, vcovmatrix_holder);

    gsl_matrix_sub(Ktest, vcovmatrix_holder);

    int k=0;
    int l=0;
    int vcov_counter=0;

    for(k=0;k<yr;k++)
        for(l=0;l<yr;l++)
    {
        vcov_matrix_out[vcov_counter++] = gsl_matrix_get(Ktest, k, l);
    }

}

/*
Perform matrix inversion for a given square matrix

@A = matrix to be inverted
@A_inverse = blank matrix to write inversion to
@N = nrows of matrix

*/

int matrix_inverse(gsl_matrix *A, gsl_matrix *A_inverse, int N)
{

    int signum = 0;
    gsl_permutation *p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_invert(A, p, A_inverse);

}

/*
Assign arbitrary vector to GSL matrix

@xrow = number of rows in new matrix
@yrow = number of columns in new vector
@K_input_matrix = data vector to be transformed to GSL matrix
@write_matrix = GSL matrix to be overwritten

*/

int matrix_assign_GSL(int *xrow, int *yrow,
           double *K_input_matrix, gsl_matrix *write_matrix)
{

    size_t nrx = xrow[0];
    size_t nry = yrow[0];
    int i, j;
    int counterin = 0;

    for (i=0;i<nrx;i++)
        for (j=0;j<nry;j++)
            {
             gsl_matrix_set(write_matrix, i, j, K_input_matrix[counterin++]);
            }
}


/*
Perform RBF kernel calculation on two matrices (x,y)


@param_in = tuning parameter for RBF Kernel
@xrow = number of rows in xmatrix
@yrow = number of rows in ymatrix
@col = number of columns in both matrices (has to be the same in both, obviously)
@input_matrix_x = x matrix (training set)
@input_matrix_y = y matrix (test set)
@matmul_matrix, @product_matrix = (xrow, xrow) blank matrix to store intermediate results during matrix operations
@kernel_matrix = (xrow, xrow) matrix that will hold output of kernel calculation
*/

int rbf_kernel(double *param_in, int *xrow, int *yrow,
               int *col, double *input_matrix_x,
               double *input_matrix_y, double *matmul_matrix,
               double *product_matrix, double *kernel_matrix)
{
    // Number of elements of mu is length(nc) and vcov is nc*nc
    size_t nrx = xrow[0];
    size_t nry = yrow[0];
    size_t nc = col[0];
    int i, j, k;
    int counterx = 0;
    int countery = 0;
    const double param = param_in[0];

    // Convert passed r matrix to gsl matrix
    gsl_matrix * x = gsl_matrix_alloc(nrx, nc);
    gsl_matrix * y = gsl_matrix_alloc(nry, nc);
    gsl_matrix * xsq = gsl_matrix_alloc(nrx, 1);
    gsl_matrix * ysq = gsl_matrix_alloc(nry, 1);
    gsl_matrix * xycombine = gsl_matrix_alloc(nrx, nry);
    gsl_matrix * mulmatrix = gsl_matrix_alloc(nrx, nry);

    // Need to assign values to matrices and calcuilate euclidean norm for each row
    // This happens for x and y matrix
    for (i=0;i<nrx;i++)
        {
        for (j=0;j<nc;j++)
            {
             gsl_matrix_set(x, i, j, input_matrix_x[counterx++]);
            }
        gsl_vector * tmp_x_vec = gsl_vector_alloc(nc);
        gsl_matrix_get_row(tmp_x_vec, x, i);
        double euclideanx = gsl_blas_dnrm2(tmp_x_vec)*gsl_blas_dnrm2(tmp_x_vec);
        gsl_matrix_set(xsq, i, 0, euclideanx);
        }

    for (i=0;i<nry;i++)
        {
        for (j=0;j<nc;j++)
            {
             gsl_matrix_set(y, i, j, input_matrix_y[countery++]);
            }
        gsl_vector * tmp_y_vec = gsl_vector_alloc(nc);
        gsl_matrix_get_row(tmp_y_vec, y, i);
        double euclideany = gsl_blas_dnrm2(tmp_y_vec)*gsl_blas_dnrm2(tmp_y_vec);
        gsl_matrix_set(ysq, i, 0, euclideany);
        }

    // Then calculate xsq + Tysq
    for (i=0;i<nrx;i++)
        for (j=0;j<nry;j++)
            {
            double tmp_entry_x = gsl_matrix_get(xsq, i, 0);
            double tmp_entry_y = gsl_matrix_get(ysq, j, 0);
            gsl_matrix_set(xycombine, i, j, tmp_entry_x+tmp_entry_y);
            }

    // Perform matrix multiplication and assign to output_matrix
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, x, y,
                  0.0, mulmatrix);

    const double scaler = 2.0;
    gsl_matrix_scale(mulmatrix, scaler);

    // Finally, assign to results_matrix for use in R
    int counterout = 0;
    int counteroutm = 0;
    for(i=0;i<nrx;i++)
        for(j=0;j<nry;j++)
        {
        matmul_matrix[counterout++] = gsl_matrix_get(mulmatrix, i, j);
        product_matrix[counteroutm++] = gsl_matrix_get(xycombine, i, j);
        }

    gsl_matrix_sub(xycombine, mulmatrix);
    const double cont = -.5*(1.0/param);
    gsl_matrix_scale(xycombine, cont);

    int counteroutcombined = 0;
    for(i=0;i<nrx;i++)
        for(j=0;j<nry;j++)
        {
        kernel_matrix[counteroutcombined++] = exp(gsl_matrix_get(xycombine, i, j));
        }
}
