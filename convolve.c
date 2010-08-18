#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <fitsio.h>

#include "util.h"
#include "memalloc.h"

//#define DEBUG 1
#ifdef DEBUG
#define DPRINTF(...) { fprintf(stderr, __VA_ARGS__); }
#else
#define DPRINTF(...) { }
#endif

int main(int argc, char *argv[])
{
  char *file1, *file2;
  char col1[4] = "1:2";
  char col2[4] = "1:2";
  char *char_columns1 = col1;
  char *char_columns2 = col2;
  int columns1[2], columns2[2], c, binary=0;
  int j,n1, N1, n2, N2, n, N;
  COORD *data1;
  COORD *data2;
  COORD *data;
  int c2 = 0;

  
  while ( (c = getopt(argc, argv, "bC:c:")) != -1) {
    switch (c) {
    case 'b':
      binary = 1;
      break;
    case 'C':
      char_columns1 = optarg;
      if (!c2) 
	char_columns2 = optarg;
      break;
    case 'c':
      char_columns2 = optarg;
      c2 = 1;
      break;
    default:
      break;
    }
  }

  if ( argc-optind == 2 ) {
    file1 = argv[optind];
    file2 = argv[optind+1];
  } else { 
    fprintf(stderr, "Usage: %s [-b] [-Cxcol:ycol] [-cxcol2:ycol2] file1 file2\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  DPRINTF("loading\n");

  /* load the data */
  intSplit(char_columns1, ":", columns1);
  intSplit(char_columns2, ":", columns2);
  n1 = loadCoords(&data1, file1, columns1[0], columns1[1]);
  n2 = loadCoords(&data2, file2, columns2[0], columns2[1]);

  /* need to pad the time-series to a power of 2 - get the power */
  int power_of_2_1;
  double power_of_2_dbl_1;
  double pow2_remainder_1;
  pow2_remainder_1 = modf( log(n1)/log(2.0), &power_of_2_dbl_1 );
  power_of_2_1 = ((int) power_of_2_dbl_1 + (pow2_remainder_1?1:0) );
  N1 = pow(2,power_of_2_1);
  int power_of_2_2;
  double power_of_2_dbl_2;
  double pow2_remainder_2;
  pow2_remainder_2 = modf( log(n2)/log(2.0), &power_of_2_dbl_2 );
  power_of_2_2 = ((int) power_of_2_dbl_2 + (pow2_remainder_2?1:0) );
  N2 = pow(2,power_of_2_2);

  if (n1>n2) {
    n = n1;
    N = N1;
    data = data1;
  } else {
    n = n2;
    N = N2;
    data = data2;
  }

  //printf("N: %ld  n: %ld\n", N, n);
  // uncomment to disable power of 2 and use mixed-radix alg. (GSL-ref S15.4)
  // - mixed algorithm still used, but gives 2^n FFT unless uncommented
  // N = n;

  /* ------------------------------------------------------ */
  /* FFT timeseries to frequency domain */
  gsl_fft_real_wavetable *Rwavetab;
  gsl_fft_halfcomplex_wavetable *HCwavetab;
  gsl_fft_real_workspace *workspace;

  /* Check the factors, do N-- if sumfactors is > 20*log2(N)  */
  /* This is *not* a traditional FFT requiring factors of 2 */
  /* -- it can use factors of 2,3,...,7 and is then N^2 for other factors */
  /* Large prime numbers are dangerous factors of N - this loop checks */
  /*  the factors and decrements N until its factors are not worse than */
  /*  about 10x the log2(N) theoretical limit for N=2^n points */
  /*   The net sacrifice in length is trivial (only a few decrements req'd) */
  /* if N=n line stays commented, it will force 2^n padding and this loop */
  /*   will succeed on the first try */
  DPRINTF("factoring\n");

  int sumfactors;
  double factor_ratio = 21.0;
  N++;
  do {

    N--;
    Rwavetab  = gsl_fft_real_wavetable_alloc(N);
    
    sumfactors = 0;
    for(j=0;j<Rwavetab->nf;j++) {
      sumfactors += (int) Rwavetab->factor[j]; 
    }
    factor_ratio = sumfactors / log2(N);
    
  } while ( factor_ratio > 20.0 );
  
  
  /* put the data into real vectors */
  double *I1, *I2, *I, *Io;
  I1 = mem_malloc(N * sizeof(double));
  I2 = mem_malloc(N * sizeof(double));
  Io = mem_malloc(N * sizeof(double));
  for (int i = 0; i<N; i++ ) { I1[i] = I2[i] = Io[i] = 0.0; }
  I = (n1>n2) ? I1 : I2;
  
  for(j=0;j<n1;j++) {  I1[j] = data1[j].y; }
  for(j=0;j<n2;j++) {  I2[j] = data2[j].y; }

  /* taper the inter one's end with cos from y[n-1] to y[0] */
  int N_taper_means = 50;
  double mean_n = 0;
  double mean_0 = 0;
  
  for(j=0;j<N_taper_means;j++) {
    mean_0 += I[j];
    mean_n += I[n-j-1];
  }
  mean_0 /= N_taper_means;
  mean_n /= N_taper_means;
  double offset = 0.5*(mean_0 + mean_n);
  double amp    = 0.5 * (mean_n - mean_0);
  double Ptaper = 2 * (N - n);
  for(j=n;j<N;j++) {  I[j] = offset + amp*cos(M_2PI*(j-n)/Ptaper); }

  /* Take the FFT */
  DPRINTF("FFT\n");

  workspace = gsl_fft_real_workspace_alloc(N);
  gsl_fft_real_transform(I1, 1, N, Rwavetab, workspace);
  gsl_fft_real_transform(I2, 1, N, Rwavetab, workspace);
  /* ------------------------------------------------------- */


  /* Take the product in the Freq dom. */
  DPRINTF("Product\n");
  Io[0] = I1[0]*I2[0];
  for (j=1;j<N/2;j++) {
    Io[2*j-1] = I1[2*j-1] * I2[2*j-1] + I1[2*j] * (I2[2*j]);
    Io[2*j]   = I1[2*j] * I2[2*j-1] - I1[2*j-1] * (I2[2*j]);
  }

  /* Inverse FFT back to time domain */
  DPRINTF("Inverting\n");
  HCwavetab = gsl_fft_halfcomplex_wavetable_alloc(N);
  gsl_fft_halfcomplex_inverse(Io, 1, N, HCwavetab, workspace);

  /* ----------------------------------------------------------- */
  /* now do the same to get the variance - not for fits output */
  if (! binary) {
    
    /* just output ascii */
    DPRINTF("output\n");
    for (j=0;j<n;j++) {  printf("%.8g  %.8g\n", data[j].x, Io[j]); }
    
  } else {
    
    /* Output - truncate the time series back to orig. length (undo padding) */
    float t_total = (float) data[n-1].x - data[0].x;
    float *fitsdata, dt = t_total/n;
    fitsfile *fptr;
    int status=0, tfields=1;
    char extname[] = "convolution";
    char *ttype[] = {"I"};
    char *tform[] = {"1E"};
    char *tunit[] = {"NA"};
    LONGLONG firstrow=1, firstelem=1;

    fitsdata = malloc(n * sizeof(float));
    for (j=0;j<n;j++) { fitsdata[j] = Io[j]; }

    char path[MAX_FILENAME];
    char basefile[MAX_FILENAME];
    strncpy(basefile, file1, MAX_FILENAME);
    basefile[strlen(basefile)-5] = '\0';
    sprintf(path, "%s.norm.fits", basefile);
    
    if (! fits_create_file(&fptr, path, &status) ) {

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
    if (status) fits_report_error(stderr, status);
    mem_free(fitsdata);

  }

  /* garbage collection */
  mem_free(I1);
  mem_free(I2);
  mem_free(Io);
  mem_free(data1);
  mem_free(data2);
  gsl_fft_real_wavetable_free(Rwavetab);
  gsl_fft_real_workspace_free(workspace);
  gsl_fft_halfcomplex_wavetable_free(HCwavetab);

  assert( show_leaks(&memlist) == 0 );
  return(EXIT_SUCCESS);
}
