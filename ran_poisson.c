#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <getopt.h>

#include "util.h"
#include "memalloc.h"

int main(int argc, char *argv[])
{
  int N;
  double mu;
  double P;

  int c, fits=0;

  while( (c  = getopt(argc, argv, "f")) != -1) {
      switch(c) {
      case 'f':
	  fits = 1;
	  break;
      default:
	  break;
      }
  }
  
  if ( argc-optind < 3 ) {
    fprintf(stderr, "Usage: %s mu N P\n", argv[0]);
    exit(EXIT_FAILURE);
  } else { 
    mu = atof(argv[optind]);
    N = atoi(argv[optind + 1]);
    P = atof(argv[optind + 2]); 
  }

  /* init the generator and get the values */
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, time(NULL));

  const double inv_N = 1.0/N;

  if (! fits) {
      /* allocate memory for values */
      unsigned int *values = mem_malloc(N * sizeof(unsigned int) );
      for (int i=0;i<N;i++) {
	  values[i] = gsl_ran_poisson(r, mu);
      }
      for (int i=0;i<N;i++) {
	  printf("%.6g %d\n", i*P*inv_N, values[i]);
      }
      mem_free(values);
  } else {
      COORD *data = mem_malloc(N*sizeof(COORD));
      for (int i = 0; i < N; ++i) {
	  data[i].x = i*P*inv_N;
	  data[i].y = gsl_ran_poisson(r, mu);
      }
      writeFitsTS("-", data, N);
      mem_free(data);
  }

  gsl_rng_free(r);
  assert( show_leaks(&memlist) == 0);

  return(EXIT_SUCCESS);
}
