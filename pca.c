#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "util.h"
#include "memalloc.h"

int main(int argc, char *argv[])
{
  int i,j;
  char *file;
  ARRAY data;
  ARRAY *datap = &data;
  int N_p;
  int N_dim;

  if ( argc < 2 ) {
    fprintf(stderr, "Usage: %s file\n", argv[0]);
    exit(EXIT_FAILURE);
  } else { 
    file = argv[1];
  }

  /* load the data */
  loadFile(file, datap);
  N_dim = datap->columns;
  if (N_dim == 0) { errQuit("No columns of data. Exiting.\n"); }


  /* normalize the data */
  double *mean     = mem_malloc(N_dim * sizeof(double));
  //double *variance = mem_malloc(N_dim * sizeof(double));

  N_p = datap->n[0];
  for (i=0; i<N_dim; i++) {
    mean[i]     = gsl_stats_mean(datap->x[i], 1, datap->n[i]);
    //variance[i] = gsl_stats_variance(datap->x[i], 1, datap->n[i]);
    
    for (j=0; j<datap->n[i]; j++) {
      //datap->x[i][j] = (datap->x[i][j] - mean[i]) / sqrt(variance[i]);
      datap->x[i][j] = (datap->x[i][j] - mean[i]);
    }

    /* truncate data in columns with excess measurements */
    if (datap->n[i] < N_p) { N_p = datap->n[i]; }

  }
  mem_free(mean);

  /* build the covariance matrix */
  double covariance;
  gsl_matrix *cov_matrix = gsl_matrix_alloc(N_dim, N_dim);
  for (i=0; i<N_dim; i++) {
    for (j=i; j<N_dim; j++) {
      covariance = gsl_stats_covariance(datap->x[i], 1, datap->x[j], 1, N_p);
      gsl_matrix_set(cov_matrix, i, j, covariance);
      gsl_matrix_set(cov_matrix, j, i, covariance);
    }
  }
  
  /* compute the eigenvectors */
  gsl_vector *eigen_vals = gsl_vector_alloc(N_dim);
  gsl_matrix *eigen_vecs = gsl_matrix_alloc(N_dim, N_dim);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N_dim);
  gsl_eigen_symmv (cov_matrix, eigen_vals, eigen_vecs, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eigen_vals, eigen_vecs, GSL_EIGEN_SORT_ABS_DESC);

  double variance_sum = 0;
  for(i=0; i<N_dim; i++) {
    variance_sum += gsl_vector_get(eigen_vals,i);
  }

  /* output */
  printf("#%9s %9s %9s   %d-Evec_components\n", 
	 "Eval", "Evalfrac", "sqrtEval", N_dim);
  double eigen_val_tmp;
  gsl_vector *x = gsl_vector_alloc(N_dim);
  gsl_vector *xprime = gsl_vector_alloc(N_dim);
  for (i=0; i<N_dim; i++) {

    eigen_val_tmp = gsl_vector_get(eigen_vals,i);
    printf("#%9.4g %9.4g %9.4g   ",  
	   eigen_val_tmp, eigen_val_tmp/variance_sum, sqrt(eigen_val_tmp) );
    for (j=0; j<N_dim; j++) {
      printf("%9.4g ", gsl_matrix_get(eigen_vecs, j, i) );
    }
    printf("\n");
  }
  /* print out the transformed coordinates */
  for (i=0; i<N_p; i++) {
    
    /* put the x values in a gsl_vector */
    for (j=0; j<N_dim; j++) {
      gsl_vector_set(x, j, datap->x[j][i]);
    }
    /* do the matrix-vector product */
    gsl_blas_dgemv(CblasTrans, 1.0, eigen_vecs, x, 0, xprime);
    for (j=0; j<N_dim; j++) {
      printf("%9.4g  ", gsl_vector_get(xprime, j));
    }
    printf("\n");
  }

  gsl_vector_free(x);
  gsl_vector_free(xprime);
  gsl_vector_free(eigen_vals);
  gsl_matrix_free(eigen_vecs);

  freeArray2D(datap);
  assert( show_leaks(&memlist) == 0 );
  
  return(EXIT_SUCCESS);
}
