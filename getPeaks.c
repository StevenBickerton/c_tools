#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#include "util.h"

void usage( char *execName) {
  fprintf(stderr, "Usage: %s [-c cx:cy] [-r] infile thresh\n", execName);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {

  char *infile;
  double thresh;
  int columns[2] = { 1,2 };
  int sign = 1, c;

  while (( c = getopt(argc, argv, "rc:")) != -1 ) {
    switch (c) {
    case 'r':
      sign = -1;
      break;
    case 'c':
      intSplit(optarg, ":", columns);
      break;
    default:
      break;
    }
  }
  if ( argc - optind == 2 ) {
    infile = argv[optind];
    thresh = atof(argv[optind+1]);
  } else {
    usage(argv[0]);
  }
  columns[0]--;
  columns[1]--;
  
  int i;
  int prev_above = 0;
  double peak_value, x_peak;
  ARRAY data, *datap;
  datap = &data;

  loadFile(infile, datap);
  double *x = datap->x[columns[0]];
  double *y = datap->x[columns[1]];
  int test, test2;

  peak_value = y[0];
  x_peak = x[0];

  for (i=0; i<datap->n[columns[1]]; i++) {
    
    // if we're above the threshold
    test = (sign>0) ? (y[i] > thresh) : (y[i] < thresh);
    if (test) {
	
      //// did we just come above
      if (! prev_above) {
	peak_value = y[i];
	x_peak = x[i];
	prev_above = 1;
	
	//// were we above already
      } else {
	test2 = (sign>0) ? (y[i] > peak_value) : (y[i] < peak_value);
	if (test2) {
	  peak_value = y[i];
	  x_peak = x[i];
	}
      }
      
      // if we're below the threshold
    } else {
      
      //// did we just drop below
      if (prev_above) {
	printf("%.8g %.8g\n", x_peak, peak_value);
	prev_above = 0;
	
	//// were we below already ... do nothing
      }
    }
  }
  
// print the value in memory if we ran out of points while in a peak,
//  but were on our way down after a peak
  if ( prev_above && peak_value > y[i-1])
    printf("%.8g %.8g\n", x_peak, peak_value);
  
  exit(EXIT_SUCCESS);
}

