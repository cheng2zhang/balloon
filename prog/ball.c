#include "ball.h"



int n = 124;
double rho = 0.9;
double tp = 1.0;
double dt = 0.002;
double thdt = 0.02;
double rm = 1.5; /* cutoff for clustering */

int nequil = 10000;
int nsteps = 100000;
int nstrdf = 100;

const char *fnpos = "ball.pos";
const char *fnrdf = "rdf.dat";



int main(void)
{
  ball_t *b;
  int istep;

  b = ball_open(n, rho);
  ball_force(b);

  for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
    ball_vv(b, dt);
    ball_vrescale(b, tp, thdt);
    if ( istep <= nequil ) continue;
    if ( istep % nstrdf == 0 ) {
      ball_cnum(b, b->g, rm);
      ball_rdf(b, b->x);
    }
  }
  graph_zhist_print(b->g);
  ball_rdf_save(b, fnrdf);
  ball_writepos(b, b->x, NULL, fnpos);
  ball_close(b);
  return 0;
}
