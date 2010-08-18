#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <getopt.h>

#include "util.h"

int main(int argc, char *argv[]) {

  char *infile, c;
  COORD *data;
  int i, N;
  int stride;
  
  int verbose = 0;
  while ( (c = getopt(argc, argv, "v")) != -1 ) {
    switch(c) {
    case 'v':
      verbose = 1;
      break;
    default:
      break;
    }
  }

  if ( argc-optind != 2 ) {
    fprintf(stderr,"Usage: %s [-v] infile stride\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    infile = argv[optind];
    stride = atoi(argv[optind+1]);
  } 
  
  /* read in the datafile */
  N = loadCoords(&data, infile, 1, 2);

  int n = N / stride;
  for(i=0;i<n;i++)
    printf("%.6f  %.6f\n", data[i*stride].x, data[i*stride].y);
  
  return(EXIT_SUCCESS);
}
