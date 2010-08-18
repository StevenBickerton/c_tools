#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <getopt.h>

#include "util.h"

int main(int argc, char *argv[]) {

  char *infile, c;
  ARRAY data;
  ARRAY *datap = &data;
  int verbose = 0;
  
  while ( (c = getopt(argc, argv, "v")) != -1 ) {
    switch(c) {
    case 'v':
      verbose = 1;  // this does nothing
      break;
    default:
      break;
    }
  }

  if ( argc-optind != 1 ) {
    fprintf(stderr,"Usage: %s [-v] infile\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    infile = argv[optind];
  } 

  /* read in the datafile */
  loadFile(infile, datap);
  printArray2D(datap);

  freeArray2D(datap);
  return(EXIT_SUCCESS);
}
