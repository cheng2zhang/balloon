#ifndef SHAPER_H__
#define SHAPER_H__



/* shape controller by spherical harmonics */
typedef struct {
  int lmin;
  int lmax;
  double **c; /* coefficients */
  double **yth; /* associated legendre polynomials */
} shaper_t;



static shaper_t *shaper_open(int lmax)
{
  shaper_t *sh;
  int l, m;

  xnew(sh, 1);
  sh->lmin = 2;
  sh->lmax = lmax;
  newarr2d(sh->c, lmax + 1, 2*lmax + 2);
  newarr2d(sh->yth, lmax + 1, lmax + 1);
  for ( l = sh->lmin; l <= sh->lmax; l++ )
    for ( m = 0; m <= 2*l + 1; m++ )
      sh->c[l][m] = rand01() - 0.5;
  return sh;
}



static void shaper_close(shaper_t *sh)
{
  delarr2d(sh->c, sh->lmax + 1);
  delarr2d(sh->yth, sh->lmax + 1);
  free(sh);
}



static void shaper_getal(shaper_t *sh, double x)
{
  int l, m;

  for ( l = sh->lmin; l <= sh->lmax; l++ )
    for ( m = 0; m <= l; m++ )
      sh->yth[l][m] = Ylm0(x, l, m);
}



static double shaper_gety(shaper_t *sh, double phi)
{
  int l, m;
  double c, s, y, z = 0;

  for ( l = sh->lmin; l <= sh->lmax; l++ ) {
    for ( m = 0; m <= l; m++ ) {
      y = sh->yth[l][m];
      c = cos(phi*m);
      s = sin(phi*m);
      z += y * (sh->c[l][2*m] * c + sh->c[l][2*m+1] * s);
    }
  }
  return z;
}



static double shaper_getr(shaper_t *sh, double phi,
    double r0, double k)
{
  double r = r0 + k * shaper_gety(sh, phi);
  if ( r < 0 ) r = 0;
  return pow(r, 1./3);
}



static void shaper_lang(shaper_t *sh, double dt, double tp)
{
  int l, m;
  double y, amp;

  l = (int) ((sh->lmax + 1) * rand01());
  amp = sqrt(2 * tp * dt);
  for ( l = sh->lmin; l <= sh->lmax; l++ ) {
    for ( m = 0; m < 2 * l + 2; m++ ) {
      y = sh->c[l][m];
      /* d^2 y / d t^2 = -nabla^2 y = -l (l + 1) y */
      y += -l * (l + 1) * y * dt + randgaus() * amp;
      if ( y > 1 ) y = 1;
      else if ( y < -1 ) y = -1;
      sh->c[l][m] = y;
    }
  }
}



#endif /* SHAPER_H__ */

