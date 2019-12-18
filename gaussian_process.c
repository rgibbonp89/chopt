#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
void main()
{
    return(0)
}

void kernel()
{
    return(0)
}


int run(int *nin, double *x)
{
    int n = nin[0];
    gsl_vector * v = gsl_vector_alloc(n);

    int i;
    double a = 20;

    for (i=0;i<n;i++)
    {
      gsl_vector_set(v, i, x[i]);
    }

    gsl_vector_add_constant(v,a);

    for (i=0;i<n;i++)
    {
      x[i] = gsl_vector_get(v, i);
    } 
}

*/

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


