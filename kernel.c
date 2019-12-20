#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/*
@xrow = count of the number of rows of the input matrix x
@yrow = count of the number of rows of the input matrix y
@col = count of the number of columns of the inout matrices
@input_matrix_x = input_matrix of type double
@input_matrix_y = input_matrix of type double
@results_matrix = vector of zeros that the function will write to and that will be output
*/

int rbf_kernel(int *xrow, int *yrow, int *col, double *input_matrix_x,
double *input_matrix_y, double *results_matrix)
{
    // Number of elements of mu is length(nc) and vcov is nc*nc
    size_t nrx = xrow[0];
    size_t nry = yrow[0];
    size_t nc = col[0];
    int i, j, k;
    int counterx = 0;
    int countery = 0;
    int counterout = 0;

    // Convert passed r matrix to gsl matrix
    // Loops to assign values to the allocated matrix
    gsl_matrix * x = gsl_matrix_alloc(nrx, nc);
    gsl_matrix * y = gsl_matrix_alloc(nry, nc);
    gsl_matrix * output_matrix = gsl_matrix_alloc(nrx, nry);

    for (i=0;i<nrx;i++)
        for (j=0;j<nc;j++)
        {
             gsl_matrix_set(x, i, j, input_matrix_x[counterx++]);
        }

    for (i=0;i<nry;i++)
        for (j=0;j<nc;j++)
        {
             gsl_matrix_set(y, i, j, input_matrix_y[countery++]);
        }

    // Perform matrix multiplication and assign to output_matrix
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,
                  1.0, x, y,
                  0.0, output_matrix);

    // Finally, assign to results_matrix for use in R
    for(i=0;i<nrx;i++)
        for(j=0;j<nry;j++)
        {
        results_matrix[counterout++] = gsl_matrix_get(output_matrix, i, j);
        }
}
