/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Sun Sep 24, 2006  15:37:41 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/Canalysis  */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
/* #include <ncurses.h>  */
/* #include <unistd.h>  */
/* #include <string.h>  */

/* #define PI 3.141592654  */



void usage(char execName[]) {
    printf ("Usage: %s x\n", execName);
    exit(1);
}

int main(int argc, char *argv[])
{

  double dx = 0.0001;
  double tol = 0.0001;
  double N, f, fp;
  int max_iter = 100;

  if ( argc < 2 )
    usage(argv[0]);
  else {
    N = atof(argv[1]);
  }

  double f1 = 1.0 - 1.0/N;
  double f2 = 1.0 - 1.0/(2.0*N);

  /* one-sided */
  double x1 = 0.5;
  double x_new = x1 + 2.0*tol;
  int i = 0;

  while ( (fabs(x1-x_new) > tol)  && i<max_iter) {
    
    x1 = x_new;
    f = f1 - (1.0 - gsl_sf_erf_Q(x1));
    fp = (gsl_sf_erf_Q(x1+dx) - gsl_sf_erf_Q(x1) )/dx;    
    x_new = x1 - f/fp;

    i++;
  }
  x1 = x_new;

  i = 0;
  double x2 = x1;
  x_new = x1 + 2.0*tol;
  while ( ( fabs(x2 - x_new) > tol) && i<max_iter) {
    x2 = x_new;
    f = f2 - (1.0 - gsl_sf_erf_Q(x2));
    fp = (gsl_sf_erf_Q(x2+dx) - gsl_sf_erf_Q(x2) ) / dx;
    x_new = x2 - f/fp;
    i++;
  }

  /* output */
  printf ("# %-6s %-8s\n","1-sided", "2-sided");
  printf ("%-9.3f %-8.3f\n", x1, x2);

  return 0;

}
