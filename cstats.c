/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Wed Sep  6, 2006  17:16:45 DST */
/* Host: kuiper */
/* Working Directory: /home/bickersj/sandbox/rmsC  */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>

#include "stats.h"
#include "util.h"
#include "memalloc.h"

int main(int argc, char *argv[])
{

    // ------=-- get options and command line args ------- //
    char c, *char_columns = NULL;
    int columns[MAX_COLUMNS], n_col;
    int suppress = 0;

    // --- define an array of booleans to zero (don't print)
    //     then set them to one if option given
    //     - if no options given, set them all to one
    int prnt[NSTATS], limit_print = 0;
    for (int i=0; i<NSTATS; i++) {
      prnt[i] = 0;
    }

    while ( (c = getopt(argc, argv, "c:mMrRdsxn")) != -1) {
      switch (c) {
      case 'c':
	char_columns = optarg;
	break;
      case 'm':
	prnt[4] = 1;
	limit_print = 1;
	break;
      case 'M':
	prnt[5] = 1;
	limit_print = 1;
	break;
      case 'r':
	prnt[7] = 1;
	limit_print = 1;
	break;
      case 'R':
	prnt[8] = 1;
	limit_print = 1;
	break;
      case 'd':
	prnt[6] = 1;
	limit_print = 1;
	break;
      case 'x':
	prnt[11] = 1;
	limit_print = 1;
	break;
      case 'n':
	prnt[9] = 1;
	limit_print = 1;
	break;
      case 's':
	suppress = 1;
	break;
      default:
	break;
      }
    }

    if (limit_print == 0) {
      for (int i=0; i<NSTATS; i++) {
	prnt[i] = 1;
      }
    }


    char *file;
    if ( argc-optind == 1 ) {
      file = argv[optind];
    } else {
      fprintf(stderr, "Usage: %s [options] file\n", argv[0]);
      fprintf(stderr, "  options: -m mean, -M 3sig clipped mean\n");
      fprintf(stderr, "           -r rms,  -S 3sig clipped rms\n");
      fprintf(stderr, "           -d median\n");
      fprintf(stderr, "           -x max\n");
      fprintf(stderr, "           -n min\n");
      fprintf(stderr, "           -s suppress header/label (print num only)\n");
      exit (EXIT_FAILURE);
    }

    // ---------------   read in the data ---------------//
    ARRAY data;
    ARRAY *datap = &data;
    loadFile(file, datap);
    n_col = datap->columns;
    if (char_columns != NULL) {
      n_col = intSplit(char_columns, ":", columns);
    } else {
      for (int i=0;i<n_col;i++) { columns[i] = i+1; }
    }
    
    if (n_col==0) {
      errQuit("No columns of data.  Exiting.\n");
    }
    
    // -----------------  get the stats  ----------------//
    char *statname[NSTATS];
    STAT stats[n_col];
    getStatNames(statname);
    for (int j=0; j<n_col; j++) {
      int col = columns[j] - 1;
      stats[j] = getStats(datap->x[col], datap->n[col]);
    }

    // print column headings
    if ( ! suppress ) {
      printf ("# %6s ","");
      for (int j=0; j<n_col; j++) { printf ("%12d ", columns[j]); }
      printf ("\n");
    }

    // print the stats
    for (int i=0; i<NSTATS; i++) {
      if (prnt[i]) {
	if (! suppress) { printf ("%-8s ", statname[i]); }
	for (int j=0; j<n_col; j++) { printf ("%+12.5g ", stats[j].all[i]); }
	printf ("\n");
      }
    }

    freeArray2D(datap); // free memory for datap
    assert( show_leaks( &memlist ) == 0 );

    return 0;
}


