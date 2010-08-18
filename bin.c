/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Mon Sep 25, 2006  11:00:59 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/Canalysis  */


#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <sys/types.h>
#include <getopt.h>

#include "stats.h"
#include "util.h"
#include "memalloc.h"

void usage(char *execName) {
  printf ("Usage: %s file col binw kern(0=box,1=gau)\n", execName);
  exit(1);
}


int main(int argc, char *argv[])
{

  int N, n, i, j, l, col;
  char *datafile, c;
  int kernel;
  double resolution, binwidth;
  double low, high;
  double hwid, binlow, value, binhigh;
  int nbinlow, nbinhigh;
  double sigma;
  double negdistsqr, contrib;
  COORD *histogram;
  STAT x;
  ARRAY data;
  ARRAY *datap = &data;
  double arglow=NAN, arghigh=NAN;
  int avgcol = 0;

  while ( (c=getopt(argc, argv, "a:l:h:")) != -1) {
    switch (c) {
    case 'a':
      avgcol = atoi(optarg);
    case 'l':
      arglow = atof(optarg);
      break;
    case 'h':
      arghigh = atof(optarg);
      break;
    default:
      break;
    }
  }

  if ( argc-optind != 4 ) {
    usage(argv[0]);
  } else {
    datafile   = argv[optind];
    col        = atoi(argv[optind+1]);
    binwidth   = atof(argv[optind+2]);
    kernel     = atoi(argv[optind+3]);
  }
  col--;

  // --------------- read in the data ----------------- //
  loadFile(datafile, datap);
  N = datap->n[col];

  if ( strcmp(argv[3], "auto") == 0 ) {
    x = meanRmsClip_d(datap->x[col], N, 3.0, 0.001, 3);
    binwidth = x.rmsClip * pow( (20.0/N), 0.2);
  }

  low  = (isnan(arglow)) ? datap->min[col]  - 1.7*binwidth : arglow;
  high = (isnan(arghigh)) ? datap->max[col]  + 1.7*binwidth : arghigh;
  
  if (kernel == 0) { resolution = binwidth; }
  else             { resolution = binwidth / 10.0; }
  
  
  //
  // Calculate the dimensions of the required array and initiallize the array
  //
  
  // allocate data bin
  n = (int) ((high-low)/resolution) + 2;
  histogram = mem_malloc( n * sizeof(COORD));

  double *averages = NULL;
  if ( avgcol ) {  averages = mem_malloc( n * sizeof(double) );  }

  // initialize the histogram values
  for (i=0; i<n; i++) {
    // use bin centres for a histogram
    histogram[i].x = low + ((kernel) ? (i+0.5) : (i)) *resolution;
    histogram[i].y = 0;
    if (avgcol) {  averages[i] = 0.0; }
  }


  // Read through the input data file and determine which bins the
  // data point will fall in.  In the case of the gaussian, the amount of
  // the contribution to each relevant bin is first determined.
    
  
  // ###################################################################
  //                   Calculations for the boxcar format
  
  if ( kernel == 0 ) {
    hwid = binwidth/2.0;
    
    for (i=0; i<N; i++ ) {
      value = datap->x[col][i];
      j = (int) trunc( (value - low) / resolution );
      if (j < 0 || j >= n) { continue; }
      histogram[j].y += 1;
      if (avgcol) { averages[j] += datap->x[avgcol-1][i]; }
    }
    
  }
  
  //######################################################################
  //                   Calculations for the gaussian format
  
  sigma = binwidth/2.0;    
  if ( kernel == 1 ) {

    for (i=0; i<N; i++ ) {
      
      value = datap->x[col][i];

      if (value >= (low-3*sigma) && value <= (high+3*sigma)) {
	binlow = value - 3*sigma;
	binhigh = value + 3*sigma;
	nbinlow = (int) ((binlow - low)/resolution) + 
	  (  (binlow>low) ? (1) : (0) );
	nbinhigh = (int) ((binhigh - low)/resolution);
	
	for (l=nbinlow; l <= nbinhigh; l++) {
	  if ((l >= 0) && (l <= n)) {
	    negdistsqr = -1 * pow( (value - (low + l*resolution)), 2);
	    contrib = exp(negdistsqr/(2 * sigma*sigma));
	    histogram[l].y += contrib;
	  }
	}
      }
      
    }
    
  }



  int nstart=0, nstop=n;
  for (j=0; j<n/2; j++) {
    if (histogram[j].y>0 && nstart==0)    { nstart = j; }
    if (histogram[n-1-j].y>0 && nstop==n) { nstop = n - j; }    
    if (nstop<n && nstart>0)              { break; }
  }
  if (! isnan(arglow))   { nstart = 0; }
  if (! isnan(arghigh))  { nstop = n; }

  double density, prob;
  if (kernel == 1) {
    for (j=nstart; j<nstop; j++) {
      density = histogram[j].y / (sigma*M_SQRT2PI); 
      prob = histogram[j].y / (sigma*M_SQRT2PI*N); 
      printf ("%.6f %.6f %.6f\n", histogram[j].x+0.5*resolution, density, prob);
    }
  } else {
    for (j=nstart; j<nstop; j++) {
      printf ("%.6f %d", histogram[j].x+0.5*resolution, (int) histogram[j].y);
      if (avgcol) {
	printf (" %.6f\n", histogram[j].y ? (averages[j]/histogram[j].y) : 0);
      } else {
	printf("\n");
      }
    }
    
  }

  freeArray2D(datap);
  mem_free(histogram);

  return 0;
  
}

