#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <math.h>


int run(double *param_in, int *xrow, int *yrow, int *col,
        double *input_matrix_x, double *input_matrix_y,
        double *matmul_matrix_train, double *product_matrix_train,
        double *kernel_matrix_train, double *L_outmatrix_train)
{
    // Pass pointers themselves as arguments to rbf_kernel()
    rbf_kernel(param_in, xrow, xrow, col, input_matrix_x,
               input_matrix_x, matmul_matrix_train, product_matrix_train,
               kernel_matrix_train, L_outmatrix_train);
    L_calc(xrow, xrow, kernel_matrix_train, L_outmatrix_train);
}


/*

Calculate Cholesky decomposition and convert back to vector (rather than GSL matrix)

*/

int L_calc(int *xrow, int *yrow,
           double *L_input_matrix, double *L_output_matrix)
{

    size_t nrx = xrow[0];
    size_t nry = yrow[0];
    int i, j;
    int counterin = 0;
    int counterout = 0;
    gsl_matrix *kernel_calc = gsl_matrix_alloc(nrx, nry);

    for (i=0;i<nrx;i++)
        for (j=0;j<nry;j++)
            {
             gsl_matrix_set(kernel_calc, i, j, L_input_matrix[counterin++]);
            }

    gsl_linalg_cholesky_decomp1(kernel_calc);

    for (i=0;i<nrx;i++)
        for (j=0;j<nry;j++)
        {
        L_output_matrix[counterout++] = gsl_matrix_get(kernel_calc, i, j);
        }
}

/*
@xrow = count of the number of rows of the input matrix x
@yrow = count of the number of rows of the input matrix y
@col = count of the number of columns of the inout matrices
@input_matrix_x = input_matrix of type double
@input_matrix_y = input_matrix of type double
@results_matrix = vector of zeros that the function will write to and that will be output
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
