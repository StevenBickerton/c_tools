/*            */
/* Steven Bickerton */
/* Dept. of Astrophysical Sciences, Princeton University */
/* bick@astro.princeton.edu*/
/* Created: Fri Apr  3, 2009  15:10:43 DST */
/* Host: bender.astro.Princeton.EDU */
/* Working Directory: /Users/bick/usr/src/Canalysis  */


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include <time.h>
#include <string.h>
/* #include <math.h>  */
/* #include <fitsio.h> */
/* #include <fftw3.h> */
/* #include <gsl/gsl_errno.h> */
/* #include <gsl/gsl_rng.h> */
/* #include <gsl/gsl_randist.h> */

/* #define PI 3.141592654  */



void usage(char execName[]) {
    printf ("Usage: %s file\n", execName);
    exit(1);
}

int main(int argc, char *argv[])
{
    int n_block = 1;
    int n_lines = 36;
    int fold = 80;
    char *file;
    char c;
    int verbose = 0;

    while( (c=getopt(argc,argv,"f:b:l:v")) != -1 ) {
      switch(c) {
      case 'f':
	  fold = atoi(optarg);
	  break;
      case 'b':
	  n_block = atoi(optarg);
	  n_lines = n_block*36;
	  break;
      case 'l':
	  n_lines = atoi(optarg);
	  break;
      case 'v':
        verbose = 1;
        break;
      default:
        break;
      }
    }

    if ( argc-optind != 1 ) {
	usage(argv[0]);
    } else {
	file = argv[optind];
    }

    int block = 2880 * n_block;
    FILE *fp = fopen(file, "r");
    char *hdr = malloc(block*sizeof(char));
    fread(hdr, sizeof(char), block, fp);
    fclose(fp);

    char sub[fold];
    for (int i=0; i<n_lines; ++i) {
	snprintf(sub, sizeof(char)*fold, "%s", hdr+i*fold);
	printf("%s\n", sub);
    }
    return 0;

}
