/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Thu Oct 19, 2006  18:24:42 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/Canalysis  */


#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fitsio.h>

#include "util.h"
#include "memalloc.h"

// --------------------------------------------------------------------
void errQuit (char *msg, ...) {

  char s[MAX_LINE_LENGTH];

  va_list args;
  va_start(args,msg);
  vsprintf(s, msg, args);
  va_end(args);

  if (errno) { perror(s); }
  else       { fprintf(stderr, s); }

  exit(EXIT_FAILURE);
}

/* -------------------------------------------------------------- */
int cmp(const void *v1, const void *v2) {  
  int i;
  double diff = *(double *)v1 - *(double *)v2;

  i = (diff > 0 ) ? 1 : -1;
  if (diff == 0) { i = 0; }
  return  i;  
}

// --------------------------------------------------------------------
FILE *openfile (char *path, char *mode) {

  FILE *fp;

  if ( path == NULL || path[0] == '-'  ) {
    return stdin;
  } else if ( (fp = fopen(path, mode)) != NULL) {
    return fp;
  } else {
    errQuit("open file (%s)", path);
  }

  return NULL;  // never reached ... surpresses compiler error
}

// ----------------------------------------------------------------------
void closefile (FILE *stream) {
  if (fclose(stream) != 0) { errQuit("close file"); }
}



/* ---------------------------------------------------------------- */
/* take a string and split it into pieces at 'tok' delimiters      */
/* put the pieces into an array and return the number found       */
int split (char *s, char *tok, char *entries[]) {
  
  int i=0;

  entries[i++] = strtok(s, tok);
  while ( (entries[i]=strtok(NULL,tok)) != NULL ) { i++; }

  return i;
}


int intSplit (char *s, char *tok, int entries[]) {

  int i,n;
  char *char_entries[MAX_COLUMNS];
  n = split(s, tok, char_entries);
  for (i=0;i<n;i++) { entries[i] = atoi(char_entries[i]); }
  
  return n;
}

int doubleSplit (char *s, char *tok, double entries[]) {

  int i,n;
  char *char_entries[MAX_COLUMNS];
  n = split(s, tok, char_entries);
  for (i=0;i<n;i++) { entries[i] = atof(char_entries[i]); }
  
  return n;
}

// --------------------------------------------------------------------
// read in all lines and columns of data into an array -- allocate mem too
void loadFile(char *path, ARRAY *datap)
{
  int ncolumns, ncolumns_max, nchars, nalnum;
  char paramLine[MAX_LINE_LENGTH]; // force 1 leading white space
  char *paramtmp;
  char *space = " \t";
  double value;
  int i;
  long nlines;


  /*  handle both FITS and text */

  size_t pathlen = strlen(path);
  if ( strstr(path+pathlen-4,"fits") != NULL ) {  /* If it's a FITS files */

    fitsfile *fptr;
    int status = 0;
    int hdutype, anynul;
    LONGLONG firstrow = 1;
    LONGLONG firstelem = 1;
    float nulval = NAN;
    int j, istart=0;
    float dtf;
    
    if (! fits_open_file(&fptr, path, READONLY, &status) ) {

      if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) /* move to 2nd HDU */
	errQuit("Unable to move to HDU\n");
      
      if (hdutype != BINARY_TBL) {
	errQuit("File does not contain a binary table.\n");
      } else {
	fits_get_num_rows(fptr, &nlines, &status);
	fits_get_num_cols(fptr, &ncolumns, &status);


	/* if there's a sample time listed, fake a time column */
	if ( fits_read_key(fptr,TFLOAT, "SAMPTIME", &dtf, NULL, &status)==0 ) {

	  datap->x[0] = mem_malloc( nlines * sizeof(double) );
	  for(j=0;j<nlines;j++) {
	    datap->x[0][j]    = (double) j*dtf;
	    datap->sum[0]    += datap->x[0][j];
	    datap->sumsqr[0] += SQR(datap->x[0][j]);
	  }
	  datap->n[0]    = (int) nlines;
	  datap->ndef[0] = (int) nlines;
	  datap->min[0]  = 0;
	  datap->max[0]  = j * ((int) nlines-1);
	  istart = 1;
	} else {
	  status = 0;
	}

	datap->columns = ncolumns+istart;
	for (i=istart; i < ncolumns + istart; i++) {
	  
	  datap->x[i] = mem_malloc( nlines * sizeof(double) );
	  fits_read_col(fptr, TDOUBLE, i-istart+1, firstrow, firstelem,
			nlines, &nulval, datap->x[i], &anynul, &status);
	  
	  datap->n[i] = (int) nlines;
	  datap->ndef[i] = (int) nlines - anynul;
	  for (j=0;j<nlines;j++) {
	    if (datap->x[i][j] < datap->min[i]) {
	      datap->min[i] = datap->x[i][j];
	    }
	    if (datap->x[i][j] > datap->max[i]) {
	      datap->max[i] = datap->x[i][j];
	    }
	    datap->sum[i]    += datap->x[i][j];
	    datap->sumsqr[i] += SQR(datap->x[i][j]);
	  }
	}

      }

      fits_close_file(fptr, &status);
    }
    if (status) { fits_report_error(stderr, status); }
    
    /* If it's an ascii file */
  } else {

    // initialize some array values.
    for (i=0; i<MAX_COLUMNS; i++) {
      datap->n[i] = 0;
      datap->ndef[i] = 0;
      datap->x[i] = NULL;
      datap->sum[i] = 0.0;
      datap->sumsqr[i] = 0.0;
    }
    
    /* If it's a regular ascii data file */
    // read in the whole line
    FILE *fp = openfile(path, "r");
    nlines = ncolumns_max = 0;
    while ( fgets(paramLine, MAX_LINE_LENGTH, fp) != NULL ) {
      
      // get rid of any comment lines
      if (paramLine[0]=='#')  { continue; }
      
      // lose trailing whitespace, replace tabs with spaces
      nchars = (int) strlen(paramLine);
      nalnum = 0;
      for (i=nchars-1; i>=0; i--) {
	
	// set trailing white-space to \0
	if ( isspace(paramLine[i]) && ! nalnum) {
	  paramLine[i] = '\0';
	} else {
	  nalnum++;
	}

      }
      
      // if the whole line was space skip it (ie. lose blank lines)
      if ( ! nalnum ) {	continue; }
      
      // split on white space
      ncolumns = 0;
      while ( (paramtmp=strtok( (ncolumns?NULL:paramLine), space)) != NULL ) {
	
	// get more memory for this column and initialize to NAN
	if (datap->x[ncolumns] == NULL) {
	  datap->x[ncolumns] = mem_malloc( (nlines + 1) * sizeof(double));
	} else {
	  datap->x[ncolumns] = mem_realloc(datap->x[ncolumns], 
					   (nlines + 1) * sizeof(double));
	}
	datap->x[ncolumns][nlines] = NAN;  // bad data flag
	datap->n[ncolumns] = nlines+1;     // new array size
	
	// record the value, the number of valid entries and some info
	value = atof(paramtmp);
	if ( !isnan(value) && !isinf(value) ) {
	  datap->x[ncolumns][nlines] = value;
	  datap->ndef[ncolumns]++;
	  if ( value < datap->min[ncolumns] || ! nlines ) {
	    datap->min[ncolumns] = value;
	  }
	  if (value > datap->max[ncolumns]  || ! nlines ) {
	    datap->max[ncolumns] = value;
	  }
	  datap->sum[ncolumns] += value;
	  datap->sumsqr[ncolumns] += value*value;
	}
	
	ncolumns++;
      }
      
      if (ncolumns > ncolumns_max) { ncolumns_max = ncolumns; }
      
      nlines++;
    }
    closefile(fp); 

    // if ncolumns_max is still zero, we got nothin ...
    if ( ! ncolumns_max ) {
      fprintf(stderr, "Warning:  No data found in file: %s\n", path);
    }
    
    datap->columns = ncolumns_max;
    
  }

}



int loadCoords (COORD **coord, char *filename, int xcol, int ycol) {

  ARRAY data;
  ARRAY *datap = &data;
  int i,n;

  xcol--;
  ycol--;

  loadFile(filename, datap);

  // the number of entries is shortest column
  n = datap->ndef[xcol];
  if (datap->ndef[ycol] < n) { n = datap->ndef[ycol]; }

  *coord = mem_malloc( n * sizeof(COORD) );

  for(i=0;i<n;i++) {
    (*coord)[i].x = datap->x[xcol][i];
    (*coord)[i].y = datap->x[ycol][i];
  }
  
  freeArray2D(datap);

  return n;
}


// --------------------------------------------------------------


void writeFitsTS (char *path, COORD *data, int n) {

  float t_total = (float) data[n-1].x - data[0].x;
  float *fitsdata;
  float dt = t_total/n;
  fitsfile *fptr;
  int j;

  // variables needed for cfitsio
  int status=0, tfields=1;
  char extname[] = "Flux_norm";
  char *ttype[] = {"I"};
  char *tform[] = {"1E"};
  char *tunit[] = {"NA"};
  LONGLONG firstrow=1, firstelem=1;
  
  fitsdata = mem_malloc( n * sizeof(float) );
  for (j=0;j<n;j++) { fitsdata[j] = data[j].y; }

  if (! fits_create_file(&fptr, path, &status) ) {
    
    /* write dt to the header */
    fits_create_tbl(fptr, BINARY_TBL, n, tfields, 
		    ttype, tform, tunit, extname, &status);
    
    fits_write_key(fptr, TFLOAT, "SAMPTIME", &dt, 
		   "Sampling time", &status);
    fits_write_key(fptr, TFLOAT, "DURATION", &t_total, 
		   "Duration of timeseries", &status);
    
    /* write the data */
    fits_write_col(fptr, TFLOAT, 1, firstrow, firstelem,
		   n, fitsdata, &status);
    
  }
  fits_close_file(fptr, &status);
  mem_free(fitsdata);
  if (status) { fits_report_error(stderr, status); }
  
}

void writeFitsTable (char *path, ARRAY *data, int n) {

  float *fitsdata;
  fitsfile *fptr;
  int i,j;

  // variables needed for cfitsio
  int status=0, tfields=data->columns;
  char extname[] = "EXT";
  char *tform0 = "1E";
  char *tunit0 = "NA";
  char *ttype[tfields];
  char *tform[tfields];
  char *tunit[tfields];

  for (j=0;j<tfields;j++) {
    ttype[j] = mem_malloc( 3 * sizeof(char) );
    sprintf (ttype[j], "%d",j);
    tform[j] = tform0;
    tunit[j] = tunit0;
  }
  LONGLONG firstrow=1, firstelem=1;
 
  fitsdata = mem_malloc( n * sizeof(float) );
 
  if (! fits_create_file(&fptr, path, &status) ) {
    
    /* write dt to the header */
    fits_create_tbl(fptr, BINARY_TBL, n, tfields, 
		    ttype, tform, tunit, extname, &status);
    
    fits_write_key(fptr, TSTRING, "ORIGFILE", path, 
		   "Original Filename", &status);
    
    /* write the data */
    for (j=0;j<tfields;j++) {
      for (i=0;i<n;i++) { fitsdata[i] = (float) data->x[j][i]; }
      fits_write_col(fptr, TFLOAT, j+1, firstrow, firstelem,
		     n, fitsdata, &status);
    }
    
  }
  fits_close_file(fptr, &status);
  if (status) { fits_report_error(stderr, status); }
  
}



void printArray2D (ARRAY *datap) {

  int i, j;

  for (i=0; i<datap->n[0]; i++) {
    for (j=0; j<datap->columns; j++) {
      printf("%+12.12g ", datap->x[j][i]);
    }
    printf ("\n");
  }
}

void freeArray2D (ARRAY *datap) {

  int i;
  for (i=0; i<datap->columns; i++) {
    mem_free(datap->x[i]);
  }

}

