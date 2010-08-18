#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

#include "util.h"
#include "memalloc.h"

int main(int argc, char *argv[]) {

  char *infile, c, *column_str;
  int columns[2] = {1,2};
  char outfile[MAX_FILENAME];
  COORD *datap;
  int verbose = 0;
  
  while ( (c = getopt(argc, argv, "v")) != -1 ) {
    switch(c) {
    case 'c':
      column_str = optarg;
      intSplit(column_str, ":", columns);
      break;
    case 'v':
      verbose = 1;  // this does nothing
      break;
    default:
      break;
    }
  }

  if ( argc-optind != 1 ) {
    fprintf(stderr,"Usage: %s [-v] [-c c1:c2] infile\n", argv[0]);
    exit(EXIT_FAILURE);
  } else {
    infile = argv[optind];
  }

  /* Piped input will be printed to stdout */
  if (infile[0] == '-') {
    sprintf(outfile, "-");
  } else {
    sprintf(outfile, "!%s.fits", infile);
  }

  /* read in the datafile */
  int n = loadCoords(&datap, infile, columns[0], columns[1]);
  writeFitsTS(outfile, datap, n);

  mem_free(datap);
  
  assert( show_leaks(&memlist) == 0);
  return(EXIT_SUCCESS);
}
