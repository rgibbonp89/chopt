#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

int run(int *nin_r, int *nin_c, double *vec, double *mat)
{
    //  Initialize dimensions and allocate memory for vector and matrix
    int nr = nin_r[0];
    int nc = nin_c[0];
    gsl_matrix * m = gsl_matrix_alloc(nr, nc);
    gsl_vector * v = gsl_vector_alloc(nr);

    // Create indices and counters
    int i, j, counter;
    double a = 20;
    int c_index = 0;

    // Loops to assign values to the allocated matrix
    for (i=0;i<nr;i++)
    {
        for (j=0;j<nc;j++)
        {
             gsl_matrix_set(m, i, j, mat[counter++]);
        }
    }

    // As our main operation, add a constant to all elements of the matrix
    //   and then take the first column and assign it to a column view
    gsl_matrix_add_constant(m, a);
    gsl_vector_view column = gsl_matrix_column(m, c_index);

    for (i=0;i<nr;i++)
    {
    // The only object that gets changed is the vector we fed in as a pointer
    // as everything else is internal to the C function
      vec[i] = gsl_vector_get(&column.vector, i);
    }
}


int kernel(int *row, int *col, double *input_matrix, double *results_vector)
{
   /*
    @row = count of the number of rows of the input matrix
    @col = count of the number of columns of the inout matrix
    @input_matrix = input_matrix of type double
    @results_vector = vector of zeros that the function will write to and that will be output
    */
    // number of elements of mu is length(nc) and vcov is nc*nc
    size_t nr = row[0];
    size_t nc = col[0];
    int i, j, k;
    int counter = 0;

    // Vector of 0s for MVN means
    gsl_vector * mu_vector = gsl_vector_calloc(nc);
    gsl_vector_set_zero(mu_vector);

    // vcov matrix and set to identity
    gsl_matrix * vcov_matrix = gsl_matrix_calloc(nc, nc);
    gsl_matrix_set_all(vcov_matrix, 1);
    gsl_matrix_set_identity(vcov_matrix);
    FILE * f = fopen("test.dat", "wb");
    gsl_matrix_fwrite(f, vcov_matrix);
    fclose(f);

    // cholesky factor L of vcov
    gsl_linalg_cholesky_decomp1(vcov_matrix);

    // vector of results - gsl vector to store pdf output that will then be written to results_vector
    gsl_vector * res_holder = gsl_vector_alloc(nr);


    // convert passed r matrix to gsl matrix
    // Loops to assign values to the allocated matrix
    gsl_matrix * m = gsl_matrix_alloc(nr, nc);

    for (i=0;i<nr;i++)
        for (j=0;j<nc;j++)
        {
             gsl_matrix_set(m, i, j, input_matrix[counter++]);
        }
    FILE * x = fopen("mat.dat", "wb");
    gsl_matrix_fwrite(x, m);
    fclose(x);

    // Now calculate density
    for(k=0;k<nr;k++)
    {
        // vector of work (length ncol)
        gsl_vector * work_mvn = gsl_vector_calloc(nc);
        double result_tmp;
        gsl_vector_view row_tmp = gsl_matrix_row(m, k);
        gsl_ran_multivariate_gaussian_pdf (&row_tmp.vector, mu_vector, vcov_matrix, &result_tmp, work_mvn);
        results_vector[k] = result_tmp;
    }

}
