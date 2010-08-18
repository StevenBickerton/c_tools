#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<getopt.h>

#include "util.h"
#include "stats.h"


void usage(char *execName) {
  fprintf(stderr, "Usage: %s [-c x:y] [-f f_lo:f_hi] [-o oversamp] infile\n", execName);
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[]) {

  char *infile;
  int i, j, k, l;
  int columns[2] = { 1, 2 };
  int oversamp = 4;
  double freqs[2] = {0.0, 0.0}, tmp;
  int Nc = 5;
  int Nb = 5;
  int c;
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
      if (freqs[0] > freqs[1])
	errQuit("f_low must be less than f_high\n");
      break;
    case 'p':
      doubleSplit(optarg, ":", freqs);
      if (freqs[0] > freqs[1])
	errQuit("p_low must be les than p_high\n");
      tmp = freqs[0];
      freqs[0] = 1.0/freqs[1];
      freqs[1] = 1.0/tmp;
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    case 'b':
      Nb = atoi(optarg);
      break;
    case 's':
      Nc = atoi(optarg);
      break;
    default:
      break;
    }
  }
  if ( argc-optind == 1 )
    infile = argv[optind];
  else
    usage(argv[0]);
  columns[0]--;
  columns[1]--;


  /* --------------------------------------------------------------- */
  // load the data and compute mean squared value
  loadFile(infile, datap);

  double *t = datap->x[columns[0]];
  double *d = datap->x[columns[1]];
  double ti = datap->min[columns[0]];
  double tf = datap->max[columns[0]];
  double ymin = datap->min[columns[1]];
  double ymax = datap->max[columns[1]];
  double mean = datap->sum[columns[1]]/datap->n[columns[1]];
  int Nd = datap->n[columns[1]];

  int bin_number[Nd];
  int x_frac[Nd];
  double dt_min = fabs(t[1]-t[0]), tprev=t[0], variance=0.0;
  for (i=0; i<Nd; i++) {
    d[i] -= mean;
    variance += SQR(d[i]);
    x_frac[i] = (1.0-1.0e-12) * (d[i] - ymin) / (ymax - ymin);
    if ( fabs(t[i]-tprev) < dt_min && i )
      dt_min = fabs(t[i] - tprev); // WARN! ... relies on data being sequential
    tprev = t[i];
  }
  variance /= Nd - 1;


  /* ---------------------------------------------------------------  */
  // choose and initialize the frequencies
  double T = tf - ti;
  double df = 1.0/(oversamp*T);
  double f_min = (freqs[0]) ? (freqs[0]) : 1.0/(oversamp*T);
  double f_max = (freqs[1]) ? (freqs[1]) : 1.0/(2.0*dt_min);
  int Nf = (int) fabs( (f_max - f_min) / df );

  double *f;
  if (  (f = malloc(Nf*sizeof(double))) == NULL )
    errQuit("Allocation of f array\n");

  for (i=0; i<Nf; i++)
    f[i] = f_min + i*df;

  /* ---------------------------------------------------------------   */
  /* make arrays for binning of the data */
  double *sum_sqrs;
  double *mean_var;
  double *sum_data_dot_mean;
  double *p_sum_ddm;
  double *amplitude;
  double *max_variance;
  double bin_mean[Nb];
  int bin_n[Nb];

  if (  (sum_sqrs = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of sum_sqrs array\n");

  if (  (mean_var = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of mean_var array\n");

  if (  (sum_data_dot_mean = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of sum_data_dot_mean array\n");
  if (  (p_sum_ddm = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of p_sum_ddm array\n");
  if (  (amplitude = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of amplitude array\n");
  if (  (max_variance = calloc(Nf, sizeof(double))) == NULL )
    errQuit("Allocation of max_variance array\n");

  /* integrate the F distribution */
  double index_ratio = 1000, Fmax=4.0;
  int steps = (int) Fmax*index_ratio;
  double Fint[steps], dF=1.0/index_ratio, Fprev=0.0, F;
  Fint[0] = 0.0;
  for(k=1;k<steps;k++) {
    F = Fdist( (double) (k/index_ratio), Nd-Nb, Nd-1 );
    Fint[k] = Fint[k-1] + 0.5 * (F + Fprev) * dF;
    Fprev = F;
    //fprintf(stdout, "%.8f %.8g %.8g  %ld %d\n", k/index_ratio, F, Fint[k], Nd, Nb);
  }


  /* ---------------------------------------------------------------- */
  // compute

  // bin the data
  int dof, Pfrac_int;
  double P, Pfrac;
  double sum_sqrs_tmp, mean_var_tmp, model_min, model_max, model_var;
  double sum_data_dot_mean_tmp, sum_data_dot_mean_exp, max_variance_tmp;
  double integral=0;
  double sqrd_dev;

  dof = (int) Nd - Nb; 

  // loop over frequencies
  for(j=0; j<Nf; j++) {
    
    P = 1.0/f[j];

    sum_sqrs[j] = 2.0*variance;
    mean_var[j] = 0.0;
    sum_data_dot_mean[j] = 0.0;
    
    
    // loop over 'covers'
    for (l=0;l<Nc; l++) {
      
      // initilize the bin_ arrays
      for (i=0; i<Nb; i++) {
	bin_mean[i] = 0.0;
	bin_n[i] = 0;
	//red_chi2_tmp[i] = 0.0;
      }
      
      // go though and give each point a bin number // get mean for each bin
      for(i=0; i<Nd; i++) {
	Pfrac = (t[i] - ti + ((double) l / (double) Nc) * (P/Nb)) / (P / Nb);
	Pfrac_int = (int) trunc(Pfrac);
	bin_number[i] = Pfrac_int % Nb;
	
	bin_mean[bin_number[i]] += d[i];
	bin_n[bin_number[i]]++;
      }
      //   include Whittaker-Robinson variance of means.
      mean_var_tmp = 0.0;
      for (i=0; i<Nb; i++) {
	bin_mean[i] /= (bin_n[i]) ? (bin_n[i]) : 1;
	mean_var_tmp = SQR( bin_mean[i] );
      }
      mean_var_tmp /= Nb - 1;

      
      // get the sum of sqrs for the model
      sum_sqrs_tmp = 0.0;
      sum_data_dot_mean_tmp = 0.0;
      sum_data_dot_mean_exp = 0.0;
      for(i=0; i<Nd; i++) {
	sqrd_dev = SQR(d[i] - bin_mean[bin_number[i]]);
	sum_sqrs_tmp += sqrd_dev;
	sum_data_dot_mean_tmp += (d[i] * bin_mean[bin_number[i]]);
	sum_data_dot_mean_exp += SQR( bin_mean[bin_number[i]] );
	//fprintf (stderr, "%.8g %.8g  %.8g %.8g\n", d[i], bin_mean[bin_number[i]], sum_data_dot_mean_tmp, sum_data_dot_mean_exp);
      }
      sum_sqrs_tmp /= dof;
      
      if (sum_sqrs_tmp < sum_sqrs[j]) {
	sum_sqrs[j] = sum_sqrs_tmp;

	sum_data_dot_mean[j] = sum_data_dot_mean_tmp;
	p_sum_ddm[j] = (1.0/sum_data_dot_mean_exp) *
	  pow( ((Nd-1)*variance - sum_data_dot_mean[j]), (2-Nd)/2 ); 

	// get the parameters for the model
	model_min = model_max = bin_mean[0];
	model_var = max_variance[j] = SQR(bin_mean[0]);
	for (i=1; i<Nb; i++) {
	  if (bin_mean[i] < model_min)
	    model_min = bin_mean[i];
	  if (bin_mean[i] > model_max)
	    model_max = bin_mean[i];
	  
	  max_variance_tmp = SQR(bin_mean[i]);
	  model_var += max_variance_tmp;
	  if ( max_variance_tmp > max_variance[j] )
	    max_variance[j] = max_variance_tmp;
	}
	amplitude[j] = model_max - model_min;
	max_variance[j] *= Nb / model_var;
	
	//fprintf (stderr, "%.8g %.8g\n", sum_data_dot_mean_tmp, sum_data_dot_mean_exp);
	
      }
      
      if (mean_var_tmp > mean_var[j])
	mean_var[j] = mean_var_tmp;
      
    } // end cover l loop
    
    integral += p_sum_ddm[j] * df;
    
  } // end freq j loop
  
  
  /* -------------------------------------------------------- */
  //  output

  // header
  printf("#%7s %8s  ", "f", "P");
  printf (" pdm-%-6d P(pdm)-%-4d  ws-%-7d  modLombScar  ampl  max_variance", Nb, Nb, Nb);
  printf ("\n");
  
  // data
  double F_prob;
  for(j=0;j<Nf; j++) {
    
    F = sum_sqrs[j]/variance;
    F_prob = Fint[(int) trunc(F*index_ratio)];
    
    printf("%.6f %.6f   %.8g %.8g   %.8g   %.8g %.8g  %.8g %.8g\n", 
	   f[j], 1.0/f[j],   F, F_prob, 
	   1.0 - mean_var[j]/variance, 
	   sum_data_dot_mean[j], p_sum_ddm[j]/integral,
	   amplitude[j], max_variance[j]); 
  }

  free(f);
  free(sum_sqrs);
  free(mean_var);
  free(sum_data_dot_mean);
  free(p_sum_ddm);
  free(amplitude);
  free(max_variance);
  freeArray2D(datap);

  return 0;
}

