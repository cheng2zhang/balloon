#include "ball.h"

int n = 32;
double rho = 0.7;
double tp = 1.0;
const char *fnpos = "ball.pos";
int nequil = 1000;
int nsteps = 10000;
double dt = 0.002;
double thdt = 0.02;


int main(void)
{
  ball_t *b;
  int istep;

  b = ball_open(n, rho);
  b->epot = ball_force(b, b->x, b->f, &b->vir);

  for ( istep = 1; istep <= nequil + nsteps; istep++ ) {
    ball_vv(b, dt);
    ball_vrescale(b, tp, thdt);
    if ( istep <= nequil ) continue;
  }
  ball_writepos(b, b->x, NULL, fnpos);
  ball_close(b);
  return 0;
}
