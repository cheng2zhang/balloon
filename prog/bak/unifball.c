
/* splot "ball.pos" w p pt 7 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define PI 3.1415926535897932846



/* list perfect uniform-ball numbers */
static void ubnum(void)
{
  int i, nth, n, nb;
  double th;

  for ( nth = 1; nth <= 20; nth++ ) {
    n = nb = 0;
    for ( i = 0; i < nth; i++ ) {
      th = PI * (i + .5) / nth;
      n += (int) (2 * nth * sin(th));
      nb += (int) (2 * nth * sin(th) + .5);
    }
    printf("nth %2d, n %4d(%2d), nb %4d(%2d)\n",
        nth, n, (int) sqrt(n*PI)/2, nb, (int) sqrt(nb*PI)/2);
  }
}



/* place n points uniformly on a ball */
static void unifball(double (*x)[3], int n)
{
  int i, j, id, nth, nph;
  double y, th, ph, cth, sth, cph, sph;

  /* compute the number of divisions along theta
   * approximate solution of 2 nth / sin(PI/2/nth) = n */
  nth = (int) ( sqrt(n * PI) / 2 );

  /* we use the outer loop to avoid the case
   * that the current divisions cannot accommodate all points */
  for ( ; ; nth++ ) {
    for (id = 0, i = 0; i < nth && id < n; i++) {
      th = PI * (i + .5) / nth;
      cth = cos(th);
      sth = sin(th);
      /* determine the number of particles along the ring */
      nph = (int) (2 * nth * sth);
      for (j = 0; j < nph && id < n; j++) {
        ph = 2 * PI * (j + (i % 2) * .5) / nph;
        cph = cos(ph);
        sph = sin(ph);
        x[id][0] = sth * cph;
        x[id][1] = sth * sph;
        x[id][2] = cth;
        id++;
      }
    }
    if ( id >= n ) break;
  }
}


static int savepos(double (*x)[3], int n, const char *fn)
{
  FILE *fp;
  int i;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ )
    fprintf(fp, "%g %g %g\n", x[i][0], x[i][1], x[i][2]);
  fclose(fp);
  return 0;
}



int main(int argc, char **argv)
{
  int n = 32;
  double (*x)[3];

  if ( argc > 1 ) n = atoi(argv[1]);
  x = calloc(n, sizeof(x[0]));
  unifball(x, n);
  savepos(x, n, "ball.pos");
  ubnum();
  return 0;
}


