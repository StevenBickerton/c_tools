/*   stats.h header for libstats.c         */

#ifndef STATS_LIBRARY
#define STATS_LIBRARY

// constants
#define PI 3.141592654
#define STAT_MAX_LINE_LENGTH 80

#define NSTATS 17

#define NSIGMA 3
#define TOL 0.001
#define MAXITER 20


// macros

#define STIRLING(n) (1.0 + 1.0/(12.0*(n)) + 1.0/(288.0*(n)*(n)) + 139.0/(51840*(n)*(n)*(n)) )

#ifndef DB
#define DB(...) { fprintf(stderr, __VA_ARGS__); }
//#else
//#define DB(...) { }
#endif


/* structures */
typedef struct {
  int n;
  int ndef;
  double sum;
  double sumsqr;
  double mean;
  double meanClip;
  double median;
  double rms;
  double rmsClip;

  double min;
  double min_rms;
  int   min_i;
  double max;
  double max_rms;
  int   max_i;

  double skew;
  double skew_rms;
  double kurt;
  double kurt_rms;

  double all[NSTATS];
} STAT;


/* Functions  */

// max's and mins
int n_defined(double *values, int n);
double max_d(double *values, int n, int *index);
double min_d(double *values, int n, int *index);
int below_d (double *values, int n, double thresh, double *below, int *index, int nMax);
int above_d (double *values, int n, double thresh, double *above, int *index, int nMax);

int window (double *values, int n, double lo, double hi, double *in_range, int *index, int nMax);

// sum
double sum_d(double *values, int n);

// mean and rms  ... median, mode eventually
double mean_d(double *values, int n);
double rms_d(double *values, int n, double *sumsqr);
double median_d(double *values, int n);

double moment(double *values, int n, int r);
double moment_mean(double *values, int n, int r, double mean);
double skew_d(double *values, int n, double mean, double rms);
double kurt_d(double *values, int n, double mean, double rms);

// routines with sigma clipping
STAT meanRmsClip_d(double *values, int n, double nSigma, double tol, int iter);

double factorial (int N);
double binomCoef (int n, int k);
double binomDist (double p, int t, int n);
double poissonDist (double rate, int n);
double gamma (double n);
double Fdist (double f, int v1, int v2);
double chi2dist (double x2, int v);

void statErrQuit(char *msg, ...);

STAT getStats(double *values, int n);
void getStatNames(char *statname[NSTATS]);



#endif /* STATS_LIBRARY */
