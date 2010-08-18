/*            */
/* Steven Bickerton */
/* Dept. of Astrophysical Sciences, Princeton University */
/* bick@astro.princeton.edu*/
/* Created: Fri Oct  3, 2008  10:50:05 DST */
/* Host: bender.astro.Princeton.EDU */
/* Working Directory: /Users/bick/usr/src/Canalysis  */


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include <string.h>

#include "util.h"
#include "interpolate.h"
#include "memalloc.h"

void usage(char execName[]) {
    printf ("Usage: %s infile lo:inc:hi\n", execName);
    exit(1);
}

int main(int argc, char *argv[])
{
  char *infile, *range_str, *column_str;
  int columns[2] = {1, 2};
  char c;
  int verbose = 0;
  
  /* parse options */
  while( (c=getopt(argc,argv,"c:v")) != -1 ) {
    switch(c) {
    case 'c':
      column_str = optarg;
      intSplit(column_str, ":", columns);
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      break;
    }
  }
  
  /* parse arguments */
  if ( argc-optind != 2 ) {
    usage(argv[0]);
  } else {
    infile = argv[optind];
    range_str  = argv[optind+1];
  }
  
  /* split the range into values */
  double range[3];
  int n_inc = doubleSplit(range_str, ":", range);
  if (n_inc != 3) {
    errQuit("range must have 3 colon separated values:  lo:inc:hi\n");
  }
  if (range[0] > range[2]) {
    errQuit("'lo' must be less than 'hi' in lo:inc:hi\n");
  }
  
  /* read in the data */
  ARRAY data, *datap;
  datap = &data;
  loadFile(infile, datap);

  columns[0] -= 1;  columns[1] -= 1; /* change to zero-base */
  double *x = datap->x[columns[0]];
  double *y = datap->x[columns[1]];
  if ( datap->n[columns[0]] != datap->n[columns[1]] ) {
    errQuit("x and y columns must have the same number of values\n");
  }
  int n = datap->n[columns[0]];
  
  /* do the interpolation */
  int n_i = (range[2] - range[0]) / range[1] + 1;
  int i_start = 0;
  double *x_interp = mem_malloc( n_i * sizeof(double) );
  double *y_interp = mem_malloc( n_i * sizeof(double) );
  for (int i_i = 0; i_i < n_i; i_i++) {

    x_interp[i_i] = range[0] + i_i * range[1];

    if ( interpol(x_interp[i_i], &y_interp[i_i], x, y, n, &i_start) < 0 ) {
      fprintf(stderr, 
	      "# Warning: 'hi' not in data range, truncating at hi=%.6f\n",
	      x_interp[i_i-1]);
      n_i = i_i;
      break;
    };

  }

  /* output */
  for (int i_i=0; i_i < n_i; i_i++) {
    fprintf(stdout, "%.6f %.6f\n", x_interp[i_i], y_interp[i_i]);
  }


  /* clean up */
  mem_free(x_interp);
  mem_free(y_interp);

  return 0;
  
}
