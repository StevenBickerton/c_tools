#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<getopt.h>
#include<assert.h>
#include<float.h>

#include "util.h"
#include "memalloc.h"

#define MAX_PEAKS 100

void usage(char *execName) {
  fprintf(stderr, "Usage: %s [-c x:y] [-f f_lo:f_hi] [-o oversamp] infile\n", execName);
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[]) {

  char *infile;
  int i, j;
  int columns[2] = { 1, 2 };
  int oversamp = 4;
  int c;
  double freqs[2] = {0.0, 0.0}, tmp;
  ARRAY data, *datap;
  datap = &data;

  /* ---------------------------------------------------------------- */
  /* Command line  */
  while( (c = getopt(argc, argv, "c:f:p:o:b:s:")) != -1) {
    switch (c) {
    case 'c':
      intSplit(optarg, ":", columns);
      break;
    case 'f':
      doubleSplit(optarg, ":", freqs);
      if (freqs[0] > freqs[1]) {
	errQuit("f_low must be less than f_high\n");
      }
      break;
    case 'p':
      doubleSplit(optarg, ":", freqs);
      if (freqs[0] > freqs[1]) {
	errQuit("p_low must be les than p_high\n");
      }
      tmp = freqs[0];
      freqs[0] = 1.0/freqs[1];
      freqs[1] = 1.0/tmp;
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    default:
      break;
    }
  }
  if ( argc-optind == 1 ) {
    infile = argv[optind];
  } else {
    usage(argv[0]);
  }
  columns[0]--;
  columns[1]--;


  /* --------------------------------------------------------------- */
  // load the data and compute mean squared value
  loadFile(infile, datap);

  double *t = datap->x[columns[0]];
  double *d = datap->x[columns[1]];
  double ti = datap->min[columns[0]];
  double tf = datap->max[columns[0]];
  double mean = datap->sum[columns[1]]/datap->n[columns[1]];
  double Nd = datap->n[columns[1]];

  // get d2bar
  double Nd2bar = 0.0;
  double dt_min = fabs(t[1]-t[0]), tprev=t[0];
  for(i=0; i<Nd; i++) {
    d[i] -= mean;
    Nd2bar += d[i]*d[i];
    if ( fabs(t[i]-tprev) < dt_min && i ) {
      dt_min = fabs(t[i] - tprev); // WARN! ... relies on data being sequential
    }
    tprev = t[i];
  }

  
  /* ---------------------------------------------------------------  */
  // choose and initialize the frequencies
  double T = tf - ti;
  double df = 1.0/(oversamp*T);
  double f_min = (freqs[0]) ? (freqs[0]) : 1.0/(oversamp*T);
  double f_max = (freqs[1]) ? (freqs[1]) : 1.0/(2.0*dt_min);
  int Nf = (int) fabs( (f_max - f_min) / df );

  double *f, *p, *h2bar;
  f = mem_malloc( Nf * sizeof(double) );
  p = mem_malloc( Nf * sizeof(double) );
  h2bar = mem_malloc( Nf * sizeof(double) );

  for (i=0; i<Nf; i++)  { f[i] = f_min + i*df; }


  /* ---------------------------------------------------------------- */
  // compute the lombscargle periodogram
  double integral = 0.0, theta;
  double twopif, fourpif, cos2pift, sin2pift, cos4pift, sin4pift, Pdf;
  double R, I, C, S, difference;
  double Z = 1.0;

  // test each frequency
  for(j=0; j<Nf; j++) {
    
    twopif  = M_2PI * f[j];
    fourpif = M_4PI * f[j];

    // get theta
    cos4pift = sin4pift = 0.0;
    for(i=0; i<Nd; i++) {
      cos4pift += cos(fourpif*t[i]) * Z*Z;
      sin4pift += sin(fourpif*t[i]) * Z*Z;
    }
    theta = 0.5*atan2(sin4pift, cos4pift);

    // lomb-scargle p(f|D,I) at this frequency
    R = I = C = S = 0.0;
    for(i=0; i<Nd; i++) {
      cos2pift = cos(twopif*t[i] - theta);
      sin2pift = sin(twopif*t[i] - theta);
      R += d[i] * cos2pift * Z;
      I += d[i] * sin2pift * Z;
      C += cos2pift*cos2pift * Z*Z;
      S += sin2pift*sin2pift * Z*Z;
    }

    h2bar[j] = R*R/C + I*I/S;
    difference = (Nd2bar - h2bar[j]);
    p[j] = (1.0/sqrt(C*S)) * pow( difference, ((2-Nd)/2) ); 

    if ( isinf(p[j]) ) { p[j] = FLT_MAX; }

    Pdf = (j) ? p[j-1]*df + 0.5*(p[j]-p[j-1])*df : 0.0;
    integral += (j)  ?  Pdf : 0; 
    
  }
  
  /* -------------------------------------------------------- */
  //  output
  printf("# freq  period  LS-pdm  p(pdm)\n");
  for(j=0;j<Nf; j++) {
    printf("%.6f %.6f  %.8g %.8g\n", f[j], 1.0/f[j], h2bar[j], p[j]/integral);
  }

  mem_free(f);
  mem_free(p);
  mem_free(h2bar);
  freeArray2D(datap);
  assert( show_leaks(&memlist) == 0 );
  
  return 0;
}

