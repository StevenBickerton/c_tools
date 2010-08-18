/*  libstats.c library code          */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>

#include "stats.h"
#include "util.h"
#include "memalloc.h"


int n_defined(double *values, int n) {

  int i, n_defined=0;

  for(i=0; i<n; i++) {
    if (! isnan(values[i]) && ! isinf(values[i])) { n_defined++; }
  }
  return n_defined;
}

//////////////////////////////////////////////////////////////////////////
double min_d(double *values, int n, int *index) {

  int i;
  double min = values[0];

  for (i=0; i<n; i++ ) {
    if ( (! isnan(values[i]) && ( values[i]<min )) || isnan(min) ) {
      min = values[i]; 
    }
  }

  *index = i;
  return min;
}


//////////////////////////////////////////////////////////////////////////
double max_d(double *values, int n, int *index) {

  int i;
  double max = values[0];

  for (i=0; i<n; i++ ) {
    if (  (! isnan(values[i]) && values[i]>max) || isnan(max) ) {
      max = values[i];
    }
  }
  
  *index = i;
  return max;
}



//////////////////////////////////////////////////////////////////////////
int below_d (double *values, int n, double thresh, double *below, int *index, int nMax) {

  int i,found=0;
  
  for (i=0; i<n; i++) {
    if ( !isnan(values[i]) && values[i] < thresh ) {
      below[found] = values[i];
      index[found++] = i;
      if (found >= nMax) { statErrQuit("hits exceeded array size"); }
    }
  }

  return found;
}




//////////////////////////////////////////////////////////////////////////
int above_d (double *values, int n, double thresh, double *above, int *index, int nMax) {

  int i,found=0;
  
  for (i=0; i<n; i++) {
    if ( !isnan(values[i]) && values[i] > thresh ) {

      above[found] = values[i];
      index[found++] = i;
      if (found >= nMax) { statErrQuit("hits exceeded array size"); }
    }
  }

  return found;
}



///////////////////////////////////////////////////////////////////////////
int window (double *values, int n, double lo, double hi, double *in_range, int *index, int nMax) {

  int i,found=0;
  for (i=0; i<n; i++) {
    if ( !isnan(values[i]) && values[i] > lo  &&  values[i] < hi ) {
      in_range[found] = values[i];
      index[found++] = i;
      if (found >= nMax) { statErrQuit("hits exceeded array size"); }
    }
  }
  return found;
}



//////////////////////////////////////////////////////////////////////////
double mean_d(double *values, int n) {

  int i,n2=0;
  double sum=0;

  if (n==0) { return 0; }

  for (i=0; i<n ; i++) {
    if ( !isnan(values[i]) ) {
      sum += values[i];
      n2++;
    }
  }

  return (sum / n2);
}


//////////////////////////////////////////////////////////////////////////
double sum_d(double *values, int n) {

  int i;
  double sum=0;
  
  for (i=0; i<n ; i++) {
    if (!isnan(values[i]) ) { sum += values[i]; }
  }
  return (sum);
}


//////////////////////////////////////////////////////////////////////////
double rms_d(double *values, int n, double *sumsqr) {

  int i, n2;
  double sum;
  sum = 0;
  *sumsqr = 0;

  if (n==0) { return 0; }

  double mean = mean_d(values, n);
  double sum_diffsqr = 0.0;

  n2=0;
  for (i=0; i<n; i++) {
    if (!isnan(values[i])) {
      sum += values[i];
      n2++;
      sum_diffsqr += SQR(values[i] - mean);
      *sumsqr += (values[i] * values[i]);
    }
  }
  
  //double n2d = (double) n2; 		/* whoa ... burned by the n2*n2 > int */
  //return  sqrt(   *sumsqr/n2d    -   sum*sum/(n2d*n2d)   )  ;
  return sqrt( sum_diffsqr / (n2-1) );
}



//////////////////////////////////////////////////////////////////////////
double median_d(double *values, int n) {

  int i, j, nBy2 = n/2;
  double *sorted;

  if (n==0)  { return 0; }

  sorted = mem_malloc( n * sizeof(double) );
  
  j=0;
  for (i=0; i<n; i++) {
    if (!isnan(values[i])) { sorted[j++] = values[i]; }
  }
  qsort (sorted, j, sizeof(double), cmp);
       
  nBy2 = j/2;
  
  double median = 
    ( j % 2) ? sorted[nBy2] : (sorted[nBy2] + sorted[nBy2-1]) / 2.0;

  mem_free(sorted);
  return (median);

}


double moment(double *values, int n, int r) {

  int i, j;
  double *array, moment;

  if (n == 0) { return 0; }
  array = mem_malloc( n * sizeof(double) );

  j=0;
  for (i=0; i<n; i++) {
    if (!isnan(values[i])) { array[j++] = pow(values[i], r); }
  }
  moment = (j) ? (sum_d(array, j)/j) : 0;

  mem_free(array);
  return moment;
}


double moment_mean(double *values, int n, int r, double mean) {

  int i, j;
  double *array, moment_mean;

  array = mem_malloc( n * sizeof(double) );

  j=0;
  for (i=0; i<n; i++) {
    if (! isnan(values[i]) ) { array[j++] = values[i] - mean; }
  }
  moment_mean = moment(array, j, r);

  mem_free(array);
  return moment_mean;
}

double skew_d(double *values, int n, double mean, double rms) {

  double third_moment = moment_mean(values, n, 3, mean);
  double skew = (rms) ? ( third_moment / (rms*rms*rms) ) : 0;
  return skew;
}


double kurt_d (double *values, int n, double mean, double rms) {

  double fourth_moment = moment_mean(values, n, 4, mean);
  double kurt = (rms) ? ( fourth_moment / (rms*rms*rms*rms) ) : 0;
  return kurt;
}


//////////////////////////////////////////////////////////////////////////
STAT meanRmsClip_d(double *values, int n, double nSigma, double tol, int iter){

  int i, k, nGood, k_final;
  double mean, mean_last, rms;
  double *goodvalues;
  double sumsqr;
  STAT x;

  goodvalues = mem_malloc( n * sizeof(double) );

  mean = mean_d(values, n);
  rms  = rms_d (values, n, &sumsqr);
  mean_last = mean;

  k_final = 0;
  
  for(k=0; k<iter; k++) {
    
    nGood = 0;
    
    for (i=0; i< n; i++) {

      if ( ! isnan(values[i]) ) {
	if ( fabs(values[i] - mean) <= nSigma*rms ) {
	  goodvalues[nGood++] = values[i];
	}
      }
      
    }
    
    mean_last = mean;
    
    mean = mean_d(goodvalues, nGood);
    rms  = rms_d (goodvalues, nGood, &sumsqr);
    
    k_final = k;
    if (  fabs((mean_last - mean)/mean_last) < tol   ) { break; }
    
  }
  
  x.meanClip = mean;
  x.rmsClip  = rms;
  
  mem_free(goodvalues);
  return ( x );


}



/* ------------------------------------------------------------------- */
double factorial (int N) {

  int i;
  double Nf = (double) N;

  if (N < 0) { statErrQuit("factorial(): value must be >= zero.\n"); }

  double fact = 1.0;

  // stirling's approximation
  if (N > 50 ) {
    fact = (M_SQRT2PI * sqrt(Nf) * pow(Nf, N) * exp(-Nf)) *
      (1.0 + 1.0/(12.0*Nf) + 1.0/(288.0*Nf*Nf) + 139.0/(51840*Nf*Nf*Nf) );
    
    // if N is small, just compute it
  } else {
    for(i=N; i>0; i--) { fact *= i; }
  }

  return fact;
}


double binomCoef (int n, int k) {
  if (n<k) { statErrQuit("binomCoefs(): n must be greater than k\n"); }
  return ( (double) (factorial(n) / (factorial(k) * factorial(n-k))) ); 
}


double binomDist (double p, int t, int n) {
  
  if (n<t) { statErrQuit("binomDist(): n must be greater than t\n"); }
  if (p<0 || p>1.0) { statErrQuit("binomDist(): p must be between 0 and 1\n"); }
  
  int n_choose_t = binomCoef(n,t);
  double q = 1.0 - p;
  
  return ( (double) n_choose_t * pow(p,t) * pow(q, (n-t)));
}


double poissonDist (double rate, int n) {
  return ( pow(rate,n) * exp(-rate) / factorial(n) );
}



double gamma (double n) {
  
  int m;
  int i;
  
  double gamma;
  
  // if n is a positive integer just use factorial
  if ( n>0 && fabs(n-trunc(n))<1.0e-12  ) {
    gamma = factorial( (int) trunc(n+0.5) - 1 );
    //fprintf(stderr, "pos int %.6g %.6g\n", n, gamma);

    // if n is a negative integer
  } else if ( n<0 && fabs(n-trunc(n)) < 1.0e-12) {
    fprintf(stderr, "Warning: Gamma(x) not defined for neg. integers. Return 0\n");
    return 0;
    
    // if n is a positive half-integer
  } else if ( n>0 &&  fabs(n-0.5 - trunc(n-0.5) < 1.0e-12)  ) {
    m = (int) trunc(n);	
    gamma = M_SQRTPI / pow(2.0, m);
    for(i=1; i<=m; i++) { gamma *= 2.0 * i - 1.0; }
    
    //fprintf(stderr, "pos half-int\n");

    // if n is a negative half-integer
  } else if (n < 0  &&  fabs(n - 0.5 - trunc(n-0.5) < 1.0e-12) ) {
    m = (int) -(n - 0.5);
    gamma = M_SQRTPI * (m%2 ? (-1.0) : 1.0) * pow(2.0, m);
    for(i=1; i<=m; i++) { gamma /= 2.0 * i - 1.0; }
    //fprintf(stderr, "neg half-int\n");
    
    // if n is an arbitrary number (warning ... only good to a few digits)
  } else {
    n -= 1.0;
    gamma = M_SQRT2PI * sqrt(n) * pow(n, n) * exp(-n) *
      (1.0 + 1.0/(12.0*n) + 1.0/(288.0*n*n) + 139.0/(51840*n*n*n) );
  }

  return gamma;
}


double Fdist (double f, int n1, int n2) {

  double v1 = (double) n1;  // numerator dof
  double v2 = (double) n2;  // denominator dof

  // formula from Schaum's Mathematics handbook
  double g1 = gamma( (v1 + v2) / 2.0 );
  double g2 = gamma( v1 / 2.0 );
  double g3 = gamma( v2 / 2.0 );

  double term = pow(v1, (v1 / 2.0)) * pow(v2, (v2 / 2.0));
  double func1 = pow(f, (v1 / 2.0 - 1.0));
  double func2 = pow((v2 + v1 * f), (-(v1+v2)/2.0));
  
  double P = g1 / (g2 * g3) * term * func1 * func2;
  
  // when n1 and n2 are large, gamma() will return nan (and thus P will too)
  //  this solution uses a stirling's approx for the gamma functions
  if ( isnan(P) ) {
    g1 = ( STIRLING((v1+v2)/2.0) / ( STIRLING(v1/2.0)*STIRLING(v2/2.0) )) * M_1_2SQRTPI;
    g2 = pow( ((v1+v2)/(v2+v1*f)), ((v1+v2)/2.0) );
    g3 = sqrt( (v1*v2)/(v1+v2) );
    P = g1 * g2 * g3 * func1;
  }
  
  //fprintf(stdout, "%.8f %.8g %.8g %g %g  %.8g %.8g %.8g  %.8g %.8g %.8g\n", 
  //	  f, P, P2, v1, v2,  g1, g2, g3, term, func1, func2);
  return  P;
}



double chi2dist (double x2, int v) {

  double term1 = 1.0/( pow(2.0, (v / 2.0)) * gamma(v / 2.0) );
  double term2 = pow(x2, ((v - 2.0) / 2.0) );
  double term3 = exp(-x2 / 2.0);
  
  return term1 * term2 * term3;
}


/* ------------------------------------------------------------------ */

void getStatNames(char *statname[NSTATS]) {

  statname[0] = "n";
  statname[1] = "ndef";
  statname[2] = "sum";
  statname[3] = "sumsqr";
  statname[4] = "mean";
  statname[5] = "meanClip";
  statname[6] = "median";
  statname[7] = "rms";
  statname[8] = "rmsClip";
  
  statname[9] = "min_val";
  statname[10] = "min_rms";
  statname[11] = "max_val";
  statname[12] = "max_rms";

  statname[13] = "skew_val";
  statname[14] = "skew_rms";
  statname[15] = "kurt_val";
  statname[16] = "kurt_rms";

}


STAT getStats(double *values, int n) {

  STAT x, stats;
  int min_i, max_i;
  double mean, rms, min, max, skew, kurt, sumsqr;

  stats.all[0]   = stats.n                 = n;
  stats.all[1]   = stats.ndef              = n_defined(values, n);
  stats.all[2]   = stats.sum               = sum_d(values, n);
  rms            = rms_d   (values, n, &sumsqr);
  stats.all[3]   = stats.sumsqr            = sumsqr;

  // mean med rms
  mean           = mean_d(values, n);
  stats.all[4]   = stats.mean              = mean;
  x              = meanRmsClip_d (values, n, 
				  NSIGMA, TOL, MAXITER);
  stats.all[5]   = stats.meanClip          = x.meanClip;
  stats.all[6]   = stats.median            = median_d(values, n);
  
  stats.all[7]   = stats.rms               = rms;
  stats.all[8]   = stats.rmsClip           = x.rmsClip;
  
  // min max
  min            = min_d(values, n, &min_i);
  stats.all[9]   = stats.min               = min;
  stats.all[10]   = stats.min_rms           = (min-mean)/rms;
  stats.min_i    = min_i;
  max            = max_d(values, n, &max_i);;
  stats.all[11]  = stats.max               = max;
  stats.all[12]  = stats.max_rms           = (max-mean)/rms;
  stats.max_i    = max_i;

  // skew and kurt
  skew           = skew_d  (values, n, mean, rms);
  stats.all[13]  = stats.skew              = skew;
  stats.all[14]  = stats.skew_rms          = skew / (sqrt(15.0/n));
  kurt           = kurt_d  (values, n, mean, rms);
  stats.all[15]  = stats.kurt              = kurt;
  stats.all[16]  = stats.kurt_rms          = (kurt-3.0) / (sqrt(96.0/n));
  
  return stats;
}



void statErrQuit (char *msg, ...) {

  char s[STAT_MAX_LINE_LENGTH];

  va_list args;
  va_start(args,msg);
  vsprintf(s, msg, args);
  va_end(args);

  if (errno) { perror(s); }
  exit(EXIT_FAILURE);
}

