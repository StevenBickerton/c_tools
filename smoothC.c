/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Mon Mar 20, 2006  16:24:30 EST */
/* Host: kuiper */
/* Working Directory: /home/bickersj/sandbox/csmooth  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "util.h"

typedef struct {
  double nsigma;
  int iter;
  double binwidth;
  int kernel;
} PARAM;

void loadParams (char *paramfile, PARAM *param);
void smooth (COORD *data, double *mean, double *rms, int n, PARAM param);

int coordcmp (const void *p1, const void *p2) {
  return (int) ( 100.f * (((COORD*)p1)->x - ((COORD*)p2)->x) );
}

int main(int argc, char *argv[])
{
  int i,n;
  double *mean, *rms, *rmsClip;
  char *infile, *char_columns;
  int columns[2];
  COORD *data;
  char paramfilename[MAX_FILENAME];
  PARAM param;

  /* the default parameter file */
  sprintf (paramfilename, "%s/.smoothparams", getenv("HOME") );
  
  if ( argc < 5 ) {
    printf ("Usage: %s infile xcol:ycol kernel[0(box)|1(gauss)] width\n", argv[0]);
    exit(1);
  } else {
    infile              = argv[1];
    char_columns        = argv[2];
    param.kernel        = atoi(argv[3]);
    param.binwidth      = atof(argv[4]);
    if ( argc >= 6 )
      sprintf (paramfilename, "%s", argv[5]);
    if (argc > 6)
      fprintf(stderr, "Extra arguments (%s ...) being ignored.\n",argv[6]);
    
  }
  
  
  loadParams(paramfilename,&param);
  
  intSplit(char_columns, ":", columns);
  n = loadCoords(&data, infile, columns[0], columns[1]);
  qsort(data, n, sizeof(COORD), coordcmp);

  // allocate the mean and rms arrays
  if ( (mean = calloc(n, sizeof(double)) ) == NULL )
    errQuit("Failed memory allocation for mean array\n");
  if ( (rms = calloc(n, sizeof(double)) ) == NULL )
    errQuit("Failed memory allocation for rms array\n");
  if ( (rmsClip = calloc(n, sizeof(double)) ) == NULL )
    errQuit("Failed memory allocation for rmsClip array\n");
  
  smooth(data,mean,rms,n,param);
  for (i=0;i<n;i++)
    rmsClip[i] = rms[i];
  
  for (i=0;i<param.iter;i++)
    smooth(data,mean,rmsClip,n,param);

  printf("# x mean rms rmsClip\n");
  for (i=0;i<n;i++) {
    printf("%.6f %.6f %.6f %.6f\n", data[i].x, mean[i], rms[i], rmsClip[i]);
  }

  return 0;
  
}


/* ***********************************************
 *
 *  load parameter file
 *
 * ********************************************** */

void loadParams (char *paramfile, PARAM *param) {

  char paramLine[MAX_LINE_LENGTH];
  FILE *fp;
  char param_name[MAX_LINE_LENGTH];
  double tmp;
  int i = 0;

  param->nsigma = 3.0;
  param->iter = 1.0;

  if ( (fp = fopen(paramfile, "r")) == NULL) {
    fprintf (stderr, "No parameter file given, using defaults\n");
  } else {
    while (fgets(paramLine, MAX_LINE_LENGTH, fp) != NULL) {

      // get rid of any comment lines or blank lines
      if ( (paramLine[0]=='#') || ( strlen(paramLine)<2 ) ) {
	continue;
      }
      
      // read the first two columns as x,y (format of DAOPHOT .coo file)
      if ( sscanf (paramLine, "%s %le", param_name, &tmp) != 2 ) {
	errQuit("reading parameter file parameter: %d",i);
      }
      
      if ( strstr(param_name, "nsigma") != NULL ) {
	param->nsigma = tmp;
      }
      if ( strstr(param_name, "iter") != NULL ) {
	param->iter = (int) tmp;
      }
      
      i++;
    }
    closefile(fp);
  }
  
}


void smooth (COORD *data, double *mean, double *rms, int n, PARAM param) {

  int sign, count;
  int i,j,k;
  double x1,y1,x2,y2, dx, halfwidth;
  double gauss;
  double root2piSigma;
  double oneOverRoot2piSigma;
  double limit;
  double tmp, tmpsqr, integral;
  int kplus, kminus;

  halfwidth = param.binwidth/2.0;
  oneOverRoot2piSigma = M_1_SQRT2PI/halfwidth;
  root2piSigma = M_SQRT2PI * halfwidth;

  // just look-up the different exp() values rather than calling exp()
  //  compute it 100 times now, rather than during each loop.
  int nexp = 100;
  double expon[nexp];
  double i2dx = 3.0*halfwidth/nexp;
  for (i=0; i<nexp; i++) {
    dx = i * i2dx;
    expon[i] = (1.0/root2piSigma) * exp(-0.5*(dx/halfwidth)*(dx/halfwidth));
  }
  
  limit = halfwidth;
  if (param.kernel==1) {
    limit = 3.0*halfwidth;
  }

  for (i=0; i<n; i++) {
      
    x1 = data[i].x;
    y1 = data[i].y;
    
    // start at the point x1,y1 and work out from there
    tmp = 0.0;
    tmpsqr = 0.0;
    integral = 0.0;
    count = 0;
    for (j=0; j<n; j++) {
      
      // k=0 we'll go left, k=1 we'll go right
      kplus = 0;
      kminus = 0;
      for (k=0; k<2; k++) {
	sign = k ? 1 : (-1);
	
	if ( ((i+sign*j) >= n) || ((i+sign*j)<0) )
	  continue;
	
	x2 = data[i+sign*j].x;
	y2 = data[i+sign*j].y;
	
	// this will keep us from counting the centre point twice;
	if ( j==0 && k==0 )
	  continue;
	
	if ( ( (y2 > (mean[i]+param.nsigma*rms[i]) ||
		y2 < (mean[i]-param.nsigma*rms[i]) ) &&
	       (mean[i] != 0.0 && rms[i] != 0.0)     )    )
	  continue;

	dx = fabs(x2-x1);
	
	
	// if points are unevenly spaced, kplus and kminus register
	//  when the outward progression has gone past _both_ ends
	if (k==1 && dx>limit)
	  kplus++;
	if (k==0 && dx>limit)
	  kminus++;
	
	
	// -----------------------------------------------------------
	// the boxcar kernel
	if (param.kernel == 0  && dx <= halfwidth ) {
	  tmp += y2;
	  tmpsqr += y2*y2;
	  count++;

	  // ----------------------------------------------------------
	  // the gaussian kernel
	
	} else if ( param.kernel == 1 && dx < 3.0*halfwidth ) {
	  
	  // too slow
	  //gauss = (1.0/root2piSigma) * 
	  // exp(-dx*dx/(2.0*halfwidth*halfwidth));
	  
	  gauss = expon[(int) (dx/i2dx)];
	  tmp += y2*gauss;
	  tmpsqr += y2*y2*gauss;
	  integral += gauss;

	}  // end if (kernel) statements
	// -----------------------------------------------------------
	
	
	
      }  // end for (k=1,0) loop 
      
      // this condition met only if we've passed outside the bin at both ends
      if (kplus && kminus)
	break;
      
      
    } // end    for (j=0..n)
    
    if (param.kernel == 0) {
      mean[i] = tmp/count;
      rms[i] = sqrt(tmpsqr/count - mean[i]*mean[i]);
    }
    if (param.kernel == 1) {
      mean[i] = tmp/integral;
      rms[i] = sqrt(tmpsqr/integral - mean[i]*mean[i]);
    }
    
  } // end      for (i=0..n)
}
