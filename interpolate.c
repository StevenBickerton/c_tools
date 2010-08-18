/*            */
/* Steven Bickerton */
/* Dept. of Astrophysical Sciences, Princeton University */
/* bick@astro.princeton.edu*/
/* Created: Fri Oct  3, 2008  11:02:08 DST */
/* Host: bender.astro.Princeton.EDU */
/* Working Directory: /Users/bick/usr/src/Canalysis  */

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
/* #include <math.h> */

/* ======================================= */
/* Simple linear interpolation */
/* ====================================== */
int interpol(double xval, double *out, double *x, double *y, int n, int *i) {
    
  /* find the highest x[i] position that is still less than xval */
  while ( x[*i+1] < xval) {
    if ( *i >= n-1 ) { return -1; }
    ++(*i); 
  }
  *out = ( (y[*i+1] - y[*i]) / (x[*i+1] - x[*i]) ) * (xval - x[*i]) + y[*i];
  
  return(0);
}

