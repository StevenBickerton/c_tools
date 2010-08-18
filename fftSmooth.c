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
  char *infile, *char_columns;
  int columns[2], c, binary=0;
  double width;
  int j,n, N;

  while ( (c = getopt(argc, argv, "b")) != -1) {
    switch (c) {
    case 'b':
      binary = 1;
      break;
    default:
      break;
    }
  }

  if ( argc-optind == 3 ) {
    infile = argv[optind];
    char_columns = argv[optind+1];
    width = atof(argv[optind+2]);
  } else { 
    fprintf(stderr, "Usage: %s [-b] infile xcol:ycol width\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  DPRINTF("loading\n");

  /* load the data */
  COORD *data;
  intSplit(char_columns, ":", columns);
  n = loadCoords(&data, infile, columns[0], columns[1]);

  /* need to pad the time-series to a power of 2 - get the power */
  int power_of_2;
  double power_of_2_dbl;
  double pow2_remainder;
  pow2_remainder = modf( log(n)/log(2.0), &power_of_2_dbl );
  power_of_2 = ((int) power_of_2_dbl + (pow2_remainder?1:0) );
  N = pow(2,power_of_2);
 
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
    for(j=0;j<Rwavetab->nf;j++) {  sumfactors += (int) Rwavetab->factor[j]; }
    factor_ratio = sumfactors / log2(N);
    
  } while ( factor_ratio > 20.0 );
  
  
  /* put the data into real vectors */
  double *I;
  I = mem_malloc(N * sizeof(double));
  for(j=0; j<n; j++) { I[j] = data[j].y; }
  for(j=n; j<N; j++) { I[j] = 0.0; }

  /* taper the end with cos from y[n-1] to y[0] */
  int N_taper_means = 50;
  double mean_n = 0;
  double mean_0 = 0;
  for(j=0;j<N_taper_means;j++) {
    mean_0 += data[j].y;
    mean_n += data[n-j-1].y;
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
  gsl_fft_real_transform(I, 1, N, Rwavetab, workspace);
  /* ------------------------------------------------------- */


  /* Take the product with the FFT of the kernel */
  DPRINTF("Product\n");
  double gauss;

  /* F( exp{-ax^2} ) =  sqrt(pi/a) * exp{-pi^2*k^2/a} */
  /*  for a = 1/(2*sigma^2), and width=2*sigma, the FT of a gaussian is:*/
  /*  1/sqrt(2*pi*sigma^2) * exp(-x^2/(2*sigma^2))   <=>  */
  /*                     exp( -(pi^2/2) * i^2 * (width/T)^2 ) */
  double gaussFreqCoeff = width/(data[n-1].x - data[0].x);
  double jd;
  gaussFreqCoeff = 0.5 * M_PI * M_PI * gaussFreqCoeff * gaussFreqCoeff;
  double jclip = 5.0*sqrt(1.0/(2.0*gaussFreqCoeff));
  for (j=1;j<N/2;j++) {
    jd = (double) j; // must cast to double or j*j=inf for N>92682
    gauss = (jd < jclip) ? exp(-jd*jd * gaussFreqCoeff) : 0.0;
    I[2*j-1]     *= gauss; // real
    I[2*j]       *= gauss; // imag
    //printf("%f\n%f\n", I[2*j], I[2*j+1]);
  }

  /* Inverse FFT back to time domain */
  DPRINTF("Inverting\n");
  HCwavetab = gsl_fft_halfcomplex_wavetable_alloc(N);
  gsl_fft_halfcomplex_inverse(I, 1, N, HCwavetab, workspace);

  
  /* ----------------------------------------------------------- */
  /* now do the same to get the variance - not for fits output */
  if (! binary) {
    double *var;
    var = mem_malloc( N * sizeof(double) );
    for(j=0; j<n; j++) { var[j] = SQR(data[j].y-I[j]); }
    for(j=n; j<N; j++) { var[j] = 0.0; } 
    
    mean_0 = mean_n = 0;
    for(j=0;j<N_taper_means;j++) {
      mean_0 += var[j];
      mean_n += var[(n-j-1)];
    }
    mean_0 /= N_taper_means;
    mean_n /= N_taper_means;
    offset = 0.5*(mean_0 + mean_n);
    amp    = 0.5 * (mean_n - mean_0);
    for(j=n;j<N;j++) {  var[j] = offset + amp*cos(M_2PI*(j-n)/Ptaper); }
    
    gsl_fft_real_transform(var, 1, N, Rwavetab, workspace);
    for (j=1;j<N/2;j++) {
      jd = (double) j; // must cast to double or j*j=inf for N>92682
      gauss = exp(-jd*jd * gaussFreqCoeff);
      var[2*j]         *= gauss; // real
      //var[2*(N-j-1)]   *= gauss; // real -ve
      var[2*j+1]       *= gauss; // cplx
      //var[2*(N-j-1)+1] *= gauss; // cplx -ve
    }
    gsl_fft_halfcomplex_inverse(var, 1, N, HCwavetab, workspace);

    /* just output ascii */
    DPRINTF("output\n");
    for (j=0;j<n;j++) {
      printf("%.8g  %.8g %.8g %.8g  %.8g\n", 
	     data[j].x, I[j], data[j].y/I[j], data[j].y, sqrt(var[j]));
    }
    
    mem_free(var);

  } else {

    char path[MAX_FILENAME];
    char basefile[MAX_FILENAME];
    strncpy(basefile, infile, MAX_FILENAME);
    basefile[strlen(basefile)-5] = '\0';
    sprintf(path, "!%s.norm.fits", basefile);

    for (j = 0; j < n; j++) {
	data[j].y = data[j].y / I[j];
    }
    
    writeFitsTS(path, data, n);

  }

  /* garbage collection */
  mem_free(I);
  mem_free(data);

  gsl_fft_real_wavetable_free(Rwavetab);
  gsl_fft_real_workspace_free(workspace);
  gsl_fft_halfcomplex_wavetable_free(HCwavetab);

  assert( show_leaks(&memlist) == 0 );

  return(EXIT_SUCCESS);
}
