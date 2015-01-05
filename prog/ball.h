#ifndef BALL_H__
#define BALL_H__



#define D 3
#include "util.h"
#include "mtrand.h"
#include "graph.h"
#include "rdf.h"



typedef struct {
  int n; /* number of particles */
  int dof; /* degrees of freedom */
  double rho; /* surface density */
  double R; /* radius */
  double area; /* surface area */

  double (*x)[D]; /* position */
  double (*v)[D]; /* velocity */
  double (*f)[D]; /* force */
  double *r2ij;
  double epot;
  double vir;
  graph_t *g;
  hist_t *rdf;
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
  xnew(b->r2ij, n * n);

  ball_inituniform(b);

  /* initialize random velocities */
  for (i = 0; i < n; i++)
    for ( d = 0; d < D; d++ )
      b->v[i][d] = randgaus();
  ball_rattle(b, b->v, b->x);

  b->g = graph_open(n);
  b->rdf = hist_open(1, 0, 2*b->R, 0.01);
  return b;
}



static void ball_close(ball_t *b)
{
  free(b->x);
  free(b->v);
  free(b->f);
  free(b->r2ij);
  graph_close(b->g);
  hist_close(b->rdf);
  free(b);
}



#define ball_energy(b) \
  b->epot = ball_energy_low(b, b->x, b->f, b->r2ij, &b->vir)

/* compute the potential energy and virial */
__inline static double ball_energy_low(ball_t *b,
    double (*x)[D], double *r2ij, double *virial)
{
  double dx[D], dr2, dr6, fs, ep, vir;
  int i, j, n = b->n;

  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      vdiff(dx, x[i], x[j]);
      dr2 = vsqr(dx);
      r2ij[i*n + j] = dr2;
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



#define ball_force(b) \
  b->epot = ball_force_low(b, b->x, b->f, b->r2ij, &b->vir)

/* compute the force and virial, return the potential energy */
__inline static double ball_force_low(ball_t *b,
    double (*x)[D], double (*f)[D], double *r2ij, double *virial)
{
  double dx[D], fi[D], dr2, dr6, fs, ep, vir;
  int i, j, n = b->n;

  for (i = 0; i < n; i++) vzero(f[i]);
  for (ep = vir = 0, i = 0; i < n - 1; i++) {
    vzero(fi);
    for (j = i + 1; j < n; j++) {
      vdiff(dx, x[i], x[j]);
      dr2 = vsqr(dx);
      r2ij[i*n + j] = dr2;
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



/* velocity verlet */
__inline static void ball_vv(ball_t *b, double dt)
{
  int i, n = b->n;
  double dth = dt * .5;

  for (i = 0; i < n; i++) { /* VV part 1 */
    vsinc(b->v[i], b->f[i], dth);
    vsinc(b->x[i], b->v[i], dt);
  }
  ball_shake(b, b->x);
  ball_force(b);
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



/* compute a graph */
__inline static void ball_mkgraph(ball_t *b, graph_t *g, double rm)
{
  int i, j, n = b->n;
  double rm2 = rm * rm;

  graph_empty(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ )
      if ( b->r2ij[i*n + j] < rm2 )
        graph_link(g, i, j);
}



/* histogram for coordination numbers */
__inline static void ball_cnum(ball_t *b, graph_t *g, double rm)
{
  ball_mkgraph(b, g, rm);
  graph_zhist(g);
}



/* accumulate the rdf */
__inline static void ball_rdf(ball_t *b, double (*x)[D])
{
  double dx[D], dr;
  int i, j, n = b->n;

  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      vdiff(dx, x[i], x[j]);
      dr = sqrt( vsqr(dx) );
      hist_add1(b->rdf, 0, dr, 1.0, HIST_VERBOSE);
    }
  }
}



/* save rdf to file */
__inline static int ball_rdf_save(const ball_t *b, const char *fn)
{
  hist_save(b->rdf, fn, HIST_ADDAHALF | HIST_KEEPLEFT | HIST_VERBOSE);
  return 0;
}



#endif /* BALL_H__ */
