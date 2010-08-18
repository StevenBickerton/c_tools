#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <getopt.h>

#include "util.h"

int main(int argc, char *argv[]) {

  char *infile, c;
  char outfile[MAX_FILENAME];
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

  /* Piped input will be printed to stdout */
  if (infile[0] == '-') {
    sprintf(outfile, "-");
  } else {
    sprintf(outfile, "%s.fits", infile);
  }

  /* read in the datafile */
  loadFile(infile, datap);
  writeFitsTable(outfile, datap, (int)datap->n[0]);

  freeArray2D(datap);
  return(EXIT_SUCCESS);
}
