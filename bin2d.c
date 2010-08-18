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
  printf ("Usage: %s [options] file c1:c2 bin1:bin2"
	  "\n"
	  "Options:\n"
	  "  -g           use gaussian smooth\n"
	  "  -x lx:hx     x-range\n"
	  "  -y ly:hy     y-range\n"
	  "  -o value     oversampling by scale \'value\'\n"
	  "\n"
	  "Example:\n"
	  "> bin2d -g file.dat 1:2 0.2:0.2    # to bin columns 1 and 2 in width=0.2 bins with gaussian\n"
	  "\n", execName);
  exit(EXIT_FAILURE);
}


void convolve(double *conv, double *hist, int nx, int ny, double xsigpix, double ysigpix) {
  
  int i, j, ik, jk;
  int nxk = (int) floor(5.0*xsigpix);
  int nyk = (int) floor(5.0*ysigpix);
  double kern[nyk*nxk];

  /* make the kernel */
  double norm = 1.0/(2.0*M_PI*xsigpix*ysigpix);
  for(j=0;j<nyk;j++) {
    for(i=0;i<nxk;i++) {
      kern[j*nxk+i] = norm * exp( -( SQR(i/xsigpix)+ SQR(j/ysigpix))/(2.0) );
    }
  }
  
  /* convolve */
  for (j=0;j<ny;j++) {
    for(i=0;i<nx;i++) {

      /* make sure there's a value to smooth */
      if (hist[j*nx+i] > 0) {

	for (jk=-nyk+1;jk<nyk;jk++) {
	  for (ik=-nxk+1;ik<nxk;ik++) {
	    if (j+jk >= 0 && j+jk < ny && i+ik >= 0 && i+ik < nx ) {
	      conv[(j+jk)*nx+(i+ik)]  += hist[j*nx+i] *kern[abs(jk)*nxk+abs(ik)];
	    }
	  }
	}

      }	/* endif */

    }
  }

}

int main(int argc, char *argv[])
{

  int N, n, n1, n2, xi, xj, i, j;
  char *datafile, c;
  double binw[2], resolution[2];
  double lo[2], hi[2], range1[2], range2[2];
  double *histogram;
  STAT xstat, ystat;
  double x, y;
  ARRAY data;
  char *char_columns = "1:2";
  char *char_binwidth = "auto";
  char *char_range1 = NULL, *char_range2 = NULL;
  int columns[2] = {0,1};
  int kern = 0;
  double oversample = 5.0;

  while ( (c=getopt(argc, argv, "gx:y:o:")) != -1) {
    switch (c) {
    case 'g':
      kern = 1;
      break;
    case 'x':
      char_range1 = optarg;
      break;
    case 'y':
      char_range2 = optarg;
      break;
    case 'o':
      oversample = atof(optarg);
      break;
    default:
      break;
    }
  }
  
  if ( argc-optind != 3 ) {
    usage(argv[0]);
  } else {
    datafile      = argv[optind];
    char_columns  = argv[optind+1];
    char_binwidth = argv[optind+2];
  }

  // --------------- read in the data ----------------- //
  intSplit(char_columns, ":", columns);
  columns[0]--;
  columns[1]--;

  loadFile(datafile, &data);
  N = data.n[columns[0]];

  /* get the binwidths */
  if ( strcmp(char_binwidth, "auto") == 0 ) {
    xstat = meanRmsClip_d(data.x[columns[0]], N, 3.0, 0.001, 3);
    binw[0] = xstat.rmsClip * pow( (20.0/N), 0.2);
    ystat = meanRmsClip_d(data.x[columns[1]], N, 3.0, 0.001, 3);
    binw[1] = ystat.rmsClip * pow( (20.0/N), 0.2);
  } else {
    doubleSplit(char_binwidth, ":", binw);
  }

  resolution[0] = (kern) ? (binw[0] / oversample) : binw[0];
  resolution[1] = (kern) ? (binw[1] / oversample) : binw[1];

  // get the low bounds if given
  if (char_range1 != NULL) {
    doubleSplit(char_range1, ":", range1);
    lo[0] = range1[0];
    hi[0] = range1[1];
  } else {
    lo[0] = data.min[columns[0]] - 1.7*binw[0];
    hi[0] = data.max[columns[0]] + 1.7*binw[0];
  }


  // get the hi bounds if given
  if (char_range2 != NULL) {
    doubleSplit(char_range2, ":", range2);
    lo[1] = range2[0];
    hi[1] = range2[1];
  } else {
    lo[1] = data.min[columns[1]] - 1.7*binw[1];
    hi[1] = data.max[columns[1]] + 1.7*binw[1];
  }
  

  // allocate data bin
  n1 = (int) ((hi[0]-lo[0])/resolution[0]) + 2;
  n2 = (int) ((hi[1]-lo[1])/resolution[1]) + 2;
  n = n1 * n2;

  histogram = mem_malloc( n * sizeof(double));

  // initialize the histogram values
  for (j=0; j<n2; j++) {
    for (i=0; i<n1; i++) {
      histogram[j*n1+i] = 0;
    }
  }




  // ###########################################################
  //                   Calculations for the boxcar format
  
  for (i=0; i<N; i++ ) {
    xi = (int) trunc( (data.x[columns[0]][i] - lo[0]) / resolution[0] );
    xj = (int) trunc( (data.x[columns[1]][i] - lo[1]) / resolution[1] );
    if (xi < 0 || xi >= n1 || xj < 0 || xj >= n2) {
      continue;
    }
    histogram[xj*n1+xi] += 1;
  }

  double *conv;
  if (kern) {
    conv = mem_malloc( n1 * n2 * sizeof(double) );
    for (i=0; i<n1*n2; i++) { conv[i] = 0.0; }
    convolve(conv, histogram, n1, n2, oversample/2.0, oversample/2.0);
  } else {
    /* needed to suppress compiler error for *conv uninitialized */
    conv = mem_malloc( 1 * sizeof(double));
  }
  
  /* output */
  double density = 0;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) {
      x = lo[0] + (i+0.5) * resolution[0];
      y = lo[1] + (j+0.5) * resolution[1];
      density = (kern) ? conv[j*n1+i] : histogram[j*n1+i];
      printf ("%.6f %.6f %g\n", x, y, density);
    }
    printf ("\n");
  }

  freeArray2D(&data);
  mem_free(histogram);
  mem_free(conv);

  return 0;
}

