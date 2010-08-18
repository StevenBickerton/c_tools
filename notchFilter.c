#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "util.h"
#include "stats.h"

#define MAX_FREQS 100


void usage(char *execName) {
  fprintf(stderr, "Usage: %s infile ff1 [ff2 ff3 ...] [cx:y] (or pp1 for periods)\n", execName);
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[]) {

  char *infile;
  int i, j, l;
  int columns[2] = { 1, 2 };
  double f[MAX_FREQS];
  int Nf = 0;
  int Nc = 5;
  int Nb = 5;
  ARRAY data, *datap;
  datap = &data;

  /* ---------------------------------------------------------------- */
  /* Command line  */
  if ( argc < 2 )
    usage(argv[0]);
  else {
    infile = argv[1];

    for(i=2; i<argc; i++) {
      
      // column numbers
      if (argv[i][0] == 'c')
	intSplit(&argv[i][1], ":", columns);
      
      // frequency range
      if (argv[i][0] == 'f')
	f[Nf++] = atof(&argv[i][1]);

      // period range
      if (argv[i][0] == 'p')
	f[Nf++] = 1.0/atof(&argv[i][1]);
      
      // Nb - the number of bins
      if (argv[i][0] == 'b')
	Nb = atoi(&argv[i][1]);
      
      // Nc - the number of covers (s for shift - c already used for 'column')
      if (argv[i][0] == 's')
	Nc = atoi(&argv[i][1]);
      
    }
    
  }
  columns[0]--;
  columns[1]--;


  /* --------------------------------------------------------------- */
  // load the data
  loadFile(infile, datap);
  double ti = data.min[columns[0]];
  double mean = data.sum[columns[1]]/data.n[columns[1]];
  int Nd = data.n[columns[1]];
  double *t = data.x[columns[0]];
  double *d = data.x[columns[1]];
  double sum_sqrs_init = data.sumsqr[columns[1]];

  // subtract off the DC and get the variance
  for(i=0; i<Nd; i++)
    d[i] -= mean;

  /* --------------------------------------------------------------- */
  /* sort the frequencies */
  //qsort (f, Nf, sizeof(double), cmp);

  /* ---------------------------------------------------------------- */
  // compute

  // bin the data
  double bin_mean[Nb][Nc];
  int bin_n[Nb][Nc];

  int Pfrac_int;
  double P, Pfrac;
  double sum_sqrs, sum_sqrs_tmp;
  int best_cover;
  int bin_number[Nd];

  // loop over frequencies
  for(j=0; j<Nf; j++) {
    
    P = 1.0/f[j];
    best_cover = 0;
    sum_sqrs = 10.0 * sum_sqrs_init;

    // initilize the bin_ arrays
    for (i=0; i<Nb; i++) {
      for (l=0; l<Nc; l++) {
	bin_mean[i][l] = 0.0;
	bin_n[i][l] = 0;
      }
    }
    
    // loop over 'covers' and get the best phase
    for (l=0;l<Nc; l++) {
      
      // go though and give each point a bin number
      for(i=0; i<Nd; i++) {
	Pfrac = (t[i] - ti + ((double) l / (double) Nc)*(P/Nb) ) / (P / Nb);
	Pfrac_int = (int) trunc(Pfrac);
	bin_number[i] = Pfrac_int % Nb;

	// get mean for each bin
	bin_mean[ bin_number[i] ][l] += d[i];
	bin_n[ bin_number[i] ][l]++;
      }
      for(i=0; i<Nb; i++)
	bin_mean[i][l] /= (bin_n[i][l]) ? (bin_n[i][l]) : 1;
      
      // get the sum of sqrs for the model
      sum_sqrs_tmp = 0.0;
      for(i=0; i<Nd; i++)
	sum_sqrs_tmp += SQR(d[i] - bin_mean[bin_number[i]][l]);
      
      if (sum_sqrs_tmp < sum_sqrs) {
	sum_sqrs = sum_sqrs_tmp;
	best_cover = l;
      }
      
    } // end cover l loop
    
    // subtract off the best cover
    for(i=0; i<Nd; i++) {
      //fprintf (stdout, "%.5f %ld %.5f %.5f %.5f %ld\n", t[i], i, d[i], bin_mean[bin_number[i]][best_cover], d[i]-bin_mean[bin_number[i]][best_cover], bin_number[i]);
      d[i] -= bin_mean[bin_number[i]][best_cover];
    }

    
  } // end freq j loop


  // add the DC component back on
  for(i=0; i<Nd; i++)
    d[i] += mean;
  
  /* -------------------------------------------------------- */
  //  output
  printArray2D(datap);
  freeArray2D(datap);
  
  return 0;
}

