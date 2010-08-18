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
  int i, N, centre, hwid;
  
  int centering = 0;
  while ( (c = getopt(argc, argv, "c")) != -1 ) {
    switch(c) {
    case 'c':
      centering = 1;
      break;
    default:
      break;
    }
  }

  if ( argc-optind != 3 ) {
    fprintf(stderr,"Usage: %s [-v] infile centre hwid\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    infile = argv[optind];
    centre = atof(argv[optind+1]);
    hwid   = atof(argv[optind+2]);
  } 
  
  /* read in the datafile */
  N = loadCoords(&data, infile, 1, 2);
  int istart = centre-hwid;
  if (istart < 0) istart = 0;
  int iend = centre+hwid;
  if (iend > N) iend = N;
  for(i=istart;i<iend;i++)
    printf("%.6f  %.6f  %d\n", 
	   (centering ? data[i].x-data[centre].x : data[i].x), 
	   data[i].y, (centering ? i-centre : i));
  
  return(EXIT_SUCCESS);
}
