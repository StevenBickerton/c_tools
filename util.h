/*            */
/* Steven Bickerton */
/* Dept. of Physics/Astronomy, McMaster University */
/* bick@physics.mcmaster.ca*/
/* Made with makeScript, Thu Oct 19, 2006  18:24:42 DST */
/* Host: kuiper */
/* Working Directory: /1/home/bickersj/sandbox/Canalysis  */


#ifndef UTILITY_LIBRARY
#define UTILITY_LIBRARY

#define MAX_FILENAME 128
#define MAX_LINE_LENGTH 16384
#define MAX_COLUMNS 2048

// these stolen from math.h
#ifndef M_E
#define M_E            2.7182818284590452354   /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E        1.4426950408889634074   /* log_2 e */
#endif

#ifndef M_LOG10E 
#define M_LOG10E       0.43429448190325182765  /* log_10 e */
#endif

#ifndef M_LN2
#define M_LN2          0.69314718055994530942  /* log_e 2 */
#endif

#ifndef M_LN10
#define M_LN10         2.30258509299404568402  /* log_e 10 */
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846  /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2         1.57079632679489661923  /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4         0.78539816339744830962  /* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI         0.31830988618379067154  /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI         0.63661977236758134308  /* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI     1.12837916709551257390  /* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */
#endif

// these not in math.h 
#ifndef M_2PI
#define M_2PI              6.2831853071795862 /* 2*pi */
#endif

#ifndef M_4PI
#define M_4PI              12.566370614359172 /* 4*pi */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI           1.7724538509055159 /* sqrt(pi) */
#endif

#ifndef M_SQRT2PI
#define M_SQRT2PI          2.5066282746310002 /* sqrt(2*pi) */
#endif

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI        0.3989422804014327 /* 1/sqrt(2*pi) */
#endif

#ifndef M_1_2SQRTPI
#define M_1_2SQRTPI        0.28209479177387814 /* 1/(2*sqrt(pi)) */
#endif



/* macros */
#define PLINE             fprintf(stderr, "%d\n", __LINE__)
#define SQR(x)            ((x)*(x))


/* structures */
typedef struct {
  double x;
  double y;
} COORD;

typedef struct {

  int n[MAX_COLUMNS];
  int ndef[MAX_COLUMNS];
  int columns;
  double *x[MAX_COLUMNS];
  double min[MAX_COLUMNS];
  double max[MAX_COLUMNS];
  double sum[MAX_COLUMNS];
  double sumsqr[MAX_COLUMNS];

} ARRAY;


/* functions */
void errQuit (char *msg, ...);

int cmp(const void *v1, const void *v2);

FILE *openfile (char *path, char *mode);
void closefile (FILE *stream) ;


int split(char *s, char *tok, char *entries[]);
int intSplit(char *s, char *tok, int entries[]);
int doubleSplit(char *s, char *tok, double entries[]);
void loadFile(char *path, ARRAY *datap);
int loadCoords(COORD **coords, char *path, int xcol, int ycol);

void writeFitsTS(char *path, COORD *data, int n);
void writeFitsTable(char *path, ARRAY *data, int n);
void printArray2D(ARRAY *datap);

void freeArray2D(ARRAY *datap);



#endif /* UTILITY_LIBRARY */
