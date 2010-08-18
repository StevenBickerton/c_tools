#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#include "util.h"
#include "memalloc.h"

int main(int argc, char *argv[]) {

  char *infile, filebase[MAX_FILENAME], outfile[MAX_FILENAME];
  int base2min, base2max, verbose=0, suppress=0;
  COORD *data;
  int N;
  char c;

  while( (c = getopt(argc, argv, "vs")) != -1) {
    switch(c) {
    case 'v':
      verbose = 1;
      break;
    case 's':
      suppress = 1;
      break;
    default:
      break;
    }
  }
  
  if ( argc-optind == 3 ) {
    infile = argv[optind];
    base2min = atof(argv[optind+1]);
    base2max = atof(argv[optind+2]);
  } else {
    fprintf(stderr,"Usage: %s [-v(erbose)] [-s(uppress remainder)] infile base2min base2max\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  strncpy(filebase, infile, MAX_FILENAME);
  filebase[ strlen(infile) - 5 ] = '\0';
  
  /* read in the datafile */
  N = loadCoords(&data, infile, 1, 2);
  
  /* output the data to files */
  int Nprinted = 0;
  int Nremaining = N;
  int base2 = (ldexp(1,base2max) > Nremaining) ? log2(Nremaining) : base2max;
  int n = ldexp(1,base2);

  int j;
  int nfile = 0;

  float t_total = (float) data[n-1].x - data[0].x;
  float *fitsdata, dt = t_total/n;
  fitsfile *fptr;
  int status=0, tfields=1;
  char extname[] = "Flux_norm";
  char *ttype[] = {"I"};
  char *tform[] = {"1E"};
  char *tunit[] = {"NA"};
  LONGLONG firstrow=1, firstelem=1;

  /* allocate temporary storage */
  fitsdata = mem_malloc( n * sizeof(float) );

  while (Nremaining) {

    if (Nremaining < ldexp(1,base2min) ) { n = Nremaining; }

    while (n > Nremaining) {
      base2--;
      n = ldexp(1, base2);
    }
    
    for (j=0;j<n;j++) {  fitsdata[j] = data[Nprinted+j].y; }
    
    /* write the fits file */
    if ( Nremaining>=ldexp(1,base2min) ||
	 ( Nremaining<ldexp(1,base2min) && ! suppress) ) {
      sprintf(outfile, "%s.%02d.fits", filebase, nfile);

      if (! fits_create_file(&fptr, outfile, &status) ) {
	
	/* write dt to the header */
	fits_create_tbl(fptr, BINARY_TBL, n, tfields, 
			ttype, tform, tunit, extname, &status);
	
	fits_write_key(fptr, TFLOAT, "SAMPTIME", &dt, 
		       "Sampling time", &status);
	fits_write_key(fptr, TFLOAT, "DURATION", &t_total, 
		       "Duration of timeseries", &status);
	
	/* write the data */
	fits_write_col(fptr, TFLOAT, 1, firstrow, firstelem,
		       n, fitsdata, &status);
	
      }
      fits_close_file(fptr, &status);
      if (status) { fits_report_error(stderr, status); }
      
      if (verbose) { fprintf (stdout, "outfile written (%d lines)\n", n); }
      
    }

    Nremaining -= n;
    Nprinted += n;
    nfile++;
  }
  
  mem_free(fitsdata);
  mem_free(data);
  assert( show_leaks(&memlist) == 0 );

  return(EXIT_SUCCESS);
}
