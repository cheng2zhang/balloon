#ifndef BALL_H__
#define BALL_H__



#define D 3
#include "util.h"
#include "mtrand.h"



typedef struct {
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  double rho; /* surface density */
  double R; /* radius */
  double area; /* surface area */

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double epot;
  double vir;
} ball_t;



/* place n points uniformly on a ball */
static void ball_inituniform(ball_t *b)
{
  int i, j, id, n = b->n, nth, nph;
  double th, ph, cth, sth, cph, sph, R = b->R;

  nth = (int) ( sqrt(PI * n) / 2 ); /* number of divisions along theta */
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
        b->x[id][0] = sth * cph * R;
        b->x[id][1] = sth * sph * R;
        b->x[id][2] = cth * R;
        id++;
      }
    }
    if ( id >= n ) break;
  }
}



/* constrain the coordinates on the surface of the ball */
static void ball_shake(ball_t *b, double (*x)[3])
{
  int i;
  double s;

  for ( i = 0; i < b->n; i++ ) {
    s = vnorm( x[i] );
    if ( s < 1e-9 ) s = 1e-9;
    vsmul(x[i], b->R / s);
  }
}



/* remove the radial components of the velocities */
static void ball_rattle(ball_t *b, double (*v)[3], double (*x)[D])
{
  int i;
  double xv, xx;

  for ( i = 0; i < b->n; i++ ) {
    xv = vdot( v[i], x[i] );
    xx = vsqr( x[i] );
    if ( xx < 1e-18 ) xx = 1e-18;
    vsinc(v[i], x[i], -xv/xx);
  }
}



/* open an LJ system */
static ball_t *ball_open(int n, double rho)
{
  ball_t *b;
  int i, d;

  xnew(b, 1);
  b->n = n;
  b->rho = rho;
  b->area = n / rho;
  b->R = sqrt( n / 4 / PI / rho );
  b->dof = n * 2;

  xnew(b->x, n);
  xnew(b->v, n);
  xnew(b->f, n);

  ball_inituniform(b);

  /* init. random velocities */
  for (i = 0; i < n; i++)
    for ( d = 0; d < D; d++ )
      b->v[i][d] = randgaus();
  ball_rattle(b, b->v, b->x);

  return b;
}



static void ball_close(ball_t *b)
{
  free(b->x);
  free(b->v);
  free(b->f);
  free(b);
}



/* compute the potential energy and virial */
__inline static double ball_energy(ball_t *b, double (*x)[D], double *virial)
{
  double dx[D], dr2, dr6, fs, ep, vir;
  int i, j, n = b->n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      vdiff(dx, x[i], x[j]);
      dr2 = vsqr(dx);
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48 * dr6 - 24); /* f.r */
      vir += fs; /* f.r */
      ep += 4 * dr6 * (dr6 - 1);
    }
  }
  if (virial) *virial = vir;
  return ep;
}



/* compute the force and virial, return the potential energy */
__inline static double ball_force(ball_t *b, double (*x)[D], double (*f)[D], double *virial)
{
  double dx[D], fi[D], dr2, dr6, fs, ep, vir;
  int i, j, n = b->n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      vdiff(dx, x[i], x[j]);
      dr2 = vsqr(dx);
      dr2 = 1 / dr2;
      dr6 = dr2 * dr2 * dr2;
      fs = dr6 * (48 * dr6 - 24); /* f.r */
      vir += fs; /* f.r */
      fs *= dr2; /* f.r / r^2 */
      vsinc(fi, dx, fs);
      vsinc(f[j], dx, -fs);
      ep += 4 * dr6 * (dr6 - 1);
    }
    vinc(f[i], fi);
  }
  if (virial) *virial = vir;
  return ep;
}



/* velocity-verlet */
__inline static void ball_vv(ball_t *b, double dt)
{
  int i, n = b->n;
  double dth = dt * .5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(b->v[i], b->f[i], dth);
    vsinc(b->x[i], b->v[i], dt);
  }
  ball_shake(b, b->x);
  b->epot = ball_force(b, b->x, b->f, &b->vir);
  for (i = 0; i < n; i++) /* VV part 2 */
    vsinc(b->v[i], b->f[i], dth);
  ball_rattle(b, b->v, b->x);
}



/* return the kinetic energy */
static double ball_ekin(double (*v)[D], int n)
{
  int i;
  double ek = 0;

  for ( i = 0; i < n; i++ )
    ek += vsqr(v[i]);
  return .5 * ek;
}



/* exact velocity rescaling thermostat */
static double ball_vrescale(ball_t *b, double tp, double dt)
{
  int i, n = b->n;
  double ek1, ek2, s, c, r, r2;

  c = (dt < 700) ? exp(-dt) : 0;
  ek1 = ball_ekin(b->v, n);
  r = randgaus();
  r2 = randchisqr(b->dof - 1);
  ek2 = ek1 + (1 - c) * ((r2 + r * r) * tp / 2 - ek1)
      + 2 * r * sqrt(c * (1 - c) * ek1 * tp / 2);
  if (ek2 < 0) ek2 = 0;
  s = sqrt(ek2/ek1);
  for (i = 0; i < n; i++) vsmul(b->v[i], s);
  return ek2;
}




/* write positions (and possibly velocities) */
__inline static int ball_writepos(ball_t *b,
    double (*x)[D], double (*v)[D], const char *fn)
{
  FILE *fp;
  int i, d;

  if (fn == NULL) fn = "ball.pos";
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %.14e\n", D, b->n, (v != NULL), b->R);
  for ( i = 0; i < b->n; i++ ) {
    for ( d = 0; d < D; d++ )
      fprintf(fp, "%.14e ", x[i][d]);
    if ( v != NULL )
      for ( d = 0; d < D; d++ )
        fprintf(fp, "%.14e ", v[i][d]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



#endif /* BALL_H__ */
