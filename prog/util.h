#ifndef UTIL_H__
#define UTIL_H__



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#if defined(Macintosh) || defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif



#define PI          3.14159265358979323846



/* Memory allocation routines */


#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %d\n", #x, (int) (n)); \
    exit(1); } }
#endif

#define newarr(x, n) { int i_; \
  xnew(x, n); \
  for ( i_ = 0; i_ < n; i_++ ) x[i_] = 0; }

#define delarr(x) free(x)

/* copy array */
#define cparr(x, y, n) { int i_; \
  for ( i_ = 0; i_ < n; i_++ ) x[i_] = y[i_]; }

/* allocate a two-dimensional array */
#define newarr2d(x, m, n) { int j_; \
  x = calloc(m, sizeof(x[0])); \
  for ( j_ = 0; j_ < m; j_++ ) { newarr( x[j_], n ); } }

/* free a two-dimensional array */
#define delarr2d(x, m) { int j_; \
  for ( j_ = 0; j_ < m; j_++ ) { delarr( x[j_] ); } \
  free(x); x = NULL; }

/* copy array */
#define cparr2d(x, y, m, n) { int j_; \
  for ( j_ = 0; j_ < m; j_++ ) cparr( x[j_], y[j_], n ); }





/* Special functions */



/* return the associated legendre polynomial P_n^m(x) */
static double alegendre(double x, int n, int m)
{
  int i, l, neg = 0;
  double y, yp, ypp, s = 1 - x*x;

  if ( m < 0 ) {
    m = -m;
    neg = 1;
  }
  if ( m > n || s < 0 ) return 0;
  s = sqrt(s);

  /* compute P_m^m(x) = (2m - 1)!! (1 - x^2)^(m/2) */
  for (yp = 1, i = 1; i <= m; i++)
    yp *= (2*i - 1)*s;
  ypp = 0;

  /* (l + 1 - m) P_{l+1}^m(x)
   *  = (2l - 1) x P_l^m(x) - (l + m) P_{l-1}^m(x)
   * start from l = m + 1 */
  for ( l = m; l < n; l++ ) {
    y = ((2*l - 1) * x * yp -  (l + m) * ypp) / (l + 1 - m);
    ypp = yp;
    yp = y;
  }

  /* handle the negative m case
   * P_n^m(x)/(n + m)! = (-)^m P_n^{-m}(x) / (n - m)! */
  if ( neg ) {
    if ( m % 2 != 0 ) yp = -yp;
    for ( i = n + m; i > n - m; i-- ) yp /= i;
  }
  return yp;
}



/* return the spherical harmonics Y_l^m(acos x, 0) */
static double Ylm0(double x, int l, int m)
{
  int i, neg = 0;
  double y, z;

  if ( m < 0 ) {
    m = -m;
    neg = 1;
  }

  /* Y_l^m(theta, phi)
   *  = (-)^m sqrt( (l - m)! / (l + m)! ) sqrt( (2 l + 1) / (4 pi) )
   *    P_l^m(cos theta) exp(i m phi)
   * Here, we'll spare the exp(i m phi) part, for phi = 0 */
  y = alegendre(x, l, m);
  if ( m % 2 != 0 ) y = -y;

  /* compute the normalization factor */
  z = 4*PI / (2*l + 1);
  for ( i = l + m; i > l - m; i-- ) z *= i;
  y /= sqrt(z);

  if ( neg ) y = -y;
  return y;
}





/* hue value (0, 1) to color */
static void huetorgb(float *c, float h)
{
  float min = 0.3f, max = 1.0f;
  float up, down, r = 0, g = 0, b = 0, f;
  int sect;

  h *= 6.f;
  sect = (int) h;
  f = (h - sect) * (max - min);
  up   = min + f;
  down = max - f;

  switch ( sect ) {
    case 0:  r = max;  g = up;   b = min;   break; /* green goes up   */
    case 1:  r = down; g = max;  b = min;   break; /* red   goes down */
    case 2:  r = min;  g = max;  b = up;    break; /* blue  goes up   */
    case 3:  r = min;  g = down; b = max;   break; /* green goes down */
    case 4:  r = up;   g = min;  b = max;   break; /* red   goes up   */
    case 5:  r = max;  g = min;  b = down;  break; /* blue  goes down */
  }
  c[0] = r; c[1] = g; c[2] = b;
}





/* String routines */



/* remove leading and trailing spaces */
__inline static char *strstrip(char *s)
{
  char *p, *q;

  /* remove trailing spaces */
  for ( p = s + strlen(s) - 1; isspace(*p); p-- )
    *p = '\0';

  /* remove leading spaces */
  for ( p = s; *p && isspace(*p); p++ )  ;
  if ( p != s )
    for ( q = s; (*q++ = *p++) != '\0'; ) ;
  return s;
}



#define strcmpfuzzy(s, t) strncmpfuzzy(s, t, INT_MAX)

/* comparison, ignoring cases, spaces and punctuations */
__inline static int strncmpfuzzy(const char *s, const char *t, int n)
{
  int is, it, i;
  const char cset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789()[]{}";

  for ( i = 0; i < n; s++, t++, i++ ) {
    while ( *s != '\0' && strchr(cset, *s) == NULL ) s++;
    while ( *t != '\0' && strchr(cset, *t) == NULL ) t++;
    is = tolower(*s);
    it = tolower(*t);
    if ( is != it ) return is - it;
    if ( *s == '\0' ) return 0;
  }
  return 0;
}



/* check if `s' starts with `t', using fuzzy comparison */
__inline static int strstartswith(const char *s, const char *t)
{
  return strncmpfuzzy(s, t, strlen(t)) == 0;
}



/* check if a string is a nonnegative integer */
__inline static int striscnum(const char *s)
{
  for ( ; *s; s++ ) if ( !isdigit(*s) ) return 0;
  return 1;
}



#endif /* UTIL_H__ */

