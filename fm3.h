#ifndef fm3_macros
#define fm3_macros

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef fm3_maxlines
#define fm3_maxlines 10000
#warning setting fm3_maxlines to default
#endif

#ifndef fm3_stream
#define fm3_stream stderr
#warning setting fm3_stream to default
#endif

int fm3_count[fm3_maxlines];
int fm3_errcode=255;
#warning setting fm3_errcode to default

int debuglevel=0;

double fm3_time;
double fm3_gettime(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)tp.tv_sec+(double)tp.tv_usec*1E-6);
}

#define DEBUG(lvl, printf_param) \
do { \
  if(debuglevel>=lvl) { \
    fprintf(fm3_stream, __FILE__ "#%d,%d: ", __LINE__, lvl); \
    fm3_printf_fm3_stream printf_param; \
    fflush(fm3_stream); \
  } \
} while(0)

int fm3_printf_fm3_stream(const char *fmt, ...)
{
  int n;
  va_list ap;

  va_start(ap, fmt);
  n=vfprintf(fm3_stream, fmt, ap);
  va_end(ap);

  return(n);
}

void fm3_msg(const char *file, int line, const char *cond, const char *date, const char *time)
{
  int i;

  fprintf(fm3_stream, "\n");
  for(i=0; i<79; i++)
    fprintf(fm3_stream, "=");

  fprintf(fm3_stream, "\n%s#%d: \"%s\"\nerror?=%s\n(%s %s)\n",
    file, line, cond, strerror(errno), date, time);
  for(i=0; i<fm3_maxlines; i++)
    if(fm3_count[i]>0) fprintf(fm3_stream, "#%d:%d  ", i, fm3_count[i]);
  fprintf(fm3_stream, "\n");
}


#define FIF(cond) \
do { \
  fm3_count[__LINE__]++; \
  if(cond) { \
    fm3_msg(__FILE__, __LINE__, #cond, __DATE__, __TIME__); \
    exit(fm3_errcode); \
  } \
} while(0)

#define FKIF(cond, printf_param) \
do { \
  fm3_count[__LINE__]++; \
  if(cond) { \
    fm3_msg(__FILE__, __LINE__, #cond, __DATE__, __TIME__); \
    fm3_printf_fm3_stream printf_param; \
    exit(fm3_errcode); \
  } \
} while(0)


void fm3_init(int errcode)
{
  int i;
  for(i=0; i<fm3_maxlines; i++)
    fm3_count[i]=0;
  fm3_errcode=errcode;
}

#define TIME(expr) \
do { \
  fm3_time=fm3_gettime(); \
  expr; \
  fprintf(fm3_stream, "%s#%d: \"%s\" %.3fs\n", __FILE__,__LINE__,#expr,fm3_gettime()-fm3_time); \
} while(0)

#endif

/*
gcc -O9 \
-Wall \
-W \
-Wtraditional \
-Wundef \
-Wshadow \
-Wbad-function-cast \
-Wcast-qual \
-Wcast-align \
-Wwrite-strings \
-Wsign-compare \
-Waggregate-return \
-Wstrict-prototypes \
-Wredundant-decls \
-Wnested-externs \
-Winline \
-Wold-style-cast \
-Woverloaded-virtual \
-Wsynth \
-Wlong-long \
$*
*/
