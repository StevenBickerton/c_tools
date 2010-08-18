/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Wed Sep  6, 2006  17:16:45 DST */
/* Host: kuiper */
/* Working Directory: /home/bickersj/sandbox/rmsC  */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>

#include "util.h"
#include "stats.h"
#include "memalloc.h"

typedef struct {
  double x;
  double y;
  double rmsLM; // LM (local mean) differences taken from local mean
  double rmsGM; // GM (global mean) diffs taken from global mean
  int v;        // degrees of freedom
} COORDRMS;


typedef struct {
  double mean;
  double rms;
} MEANRMS;


void usage(char *execName) {
    printf ("Usage: %s file width [xcol ycol]\n", execName);
    exit(1);
}

void rmsSlide (COORDRMS *data, MEANRMS *stats, double width, int n);

int main(int argc, char *argv[])
{
    int i, j, n;
    int xcol, ycol;
    double width;
    char *file;

    ARRAY data;
    ARRAY *datap = &data;
    COORDRMS *coorddata;
    MEANRMS stats;

    if (argc > 2) {
      // required args
      file = argv[1];
      width = atof(argv[2]);
      // optional args (must be two more)
      if ( argc == 5 ) {
	xcol = atoi(argv[3]) - 1;
	ycol = atoi(argv[4]) - 1;
      } else if ( argc == 4) {
	usage(argv[0]);
      } else {
	xcol = 0;
	ycol = 1;
      }
      // usage
    } else {
      usage(argv[0]);
    }

    // read in the data
    loadFile(file, datap);//xcol, ycol, data, &n);
    n = (datap->n[xcol] > datap->n[ycol]) ? datap->n[ycol] : datap->n[xcol];

    // allocate the memory
    coorddata = mem_malloc( n * sizeof(COORDRMS) );

    // transfer the data to the COORD structure
    for (j=0;j<n; j++) {
      coorddata[j].x = datap->x[xcol][j];
      coorddata[j].y = datap->x[ycol][j];
    }

    // get the global mean and rms
    double sumsqr = 0.0;
    stats.mean = mean_d(datap->x[ycol], datap->n[ycol]);
    stats.rms  = rms_d(datap->x[ycol], datap->n[ycol], &sumsqr);

    /* no inter need the original data array */
    freeArray2D(datap);

    // get the sliding std deviation
    rmsSlide(coorddata, &stats, width, n);
    
    printf ("# %8s  %10s  %10s  %10s %3s\n", 
	    "x", "rmsLM", "rmsGM", "y_orig", "v");

    for (i=0; i<n; i++) {
      printf ("%.8f  %.8f  %.8f  %.8f %3d\n", 
	      coorddata[i].x, coorddata[i].rmsLM, 
	      coorddata[i].rmsGM, coorddata[i].y, coorddata[i].v);
    }

    mem_free(coorddata);
    assert( show_leaks(&memlist) == 0 );

    return (EXIT_SUCCESS);
}


/*  Functions  */
///////////////////////////////////////////////////////////////


void rmsSlide (COORDRMS *data, MEANRMS *stats, double width, int n) {

  int i, j, i_lo, i_hi, count;
  double ysqrd, ysum, x, dysqrd, localmean;
  double halfw = width/2.0;

  for (i=0; i<n; i++) {

    /* -------------------   local mean -------------------------  */
    // get the values at point i
    x = data[i].x;
    ysum = data[i].y;
    count = 1;

    // now work our way outward from point i
    j = 1;
    i_lo = (i-j > 0) ? i-j : 0;
    i_hi = (i+j < n-1) ? i+j : n-1;

    // while within range of i, but not off the ends of the data array.
    while ( ( (i_lo>=0) && fabs(x-data[i_lo].x) < halfw )  &&
	    ( (i_hi<=n-1)  && fabs(x-data[i_hi].x) < halfw )  ) {
      
      
      // add the contribution of the point j steps below i
      if (i_lo>0 && i_lo<n) {
	ysum += data[i_lo].y;
	count++;
      } 
      // add the contribution of the point j steps above i
      if (i_hi>0 && i_hi<n) {
	ysum += data[i_hi].y;
	count++;
      }
      
      j++;
      i_lo = (i-j >= 0) ? i-j : 0;
      i_hi = (i+j <= n-1) ? i+j : n-1;

    } 
    localmean = ysum/count;

    /* ----------------------  rms  -------------------------------*/
    // get the values at point i
    x = data[i].x;
    ysqrd = (data[i].y-localmean) * (data[i].y-localmean);
    dysqrd = (data[i].y - stats->mean)*(data[i].y - stats->mean);
    count = 1;

    // now work our way outward from point i
    j = 1;
    i_lo = (i-j > 0) ? i-j : 0;
    i_hi = (i+j < n-1) ? i+j : n-1;

    // while within range of i, but not off the ends of the data array.
    while ( ( (i_lo>=0) && fabs(x-data[i_lo].x) < halfw )  &&
	    ( (i_hi<=n-1)  && fabs(x-data[i_hi].x) < halfw )  ) {
      
      
      // add the contribution of the point j steps below i
      if (i_lo>0 && i_lo<n) {
	ysqrd += (data[i_lo].y-localmean) * (data[i_lo].y-localmean);
	dysqrd   += (data[i_lo].y - stats->mean)*(data[i_lo].y - stats->mean);
	count++;
      } 
      // add the contribution of the point j steps above i
      if (i_hi>0 && i_hi<n) {
	ysqrd += (data[i_hi].y-localmean) * (data[i_hi].y-localmean);
	dysqrd   += (data[i_hi].y - stats->mean)*(data[i_hi].y - stats->mean);
	count++;
      }
      
      j++;
      i_lo = (i-j >= 0) ? i-j : 0;
      i_hi = (i+j <= n-1) ? i+j : n-1;

    } 

    data[i].rmsLM = sqrt( ysqrd/(count-1) );
    data[i].rmsGM = sqrt(  dysqrd/(count-1)  );    
    data[i].v = count-1;
  }

}
