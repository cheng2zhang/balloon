#ifndef GRAPH_H__
#define GRAPH_H__



typedef struct {
  int n;
  char *mat;
  double *zhist;
} graph_t;



static graph_t *graph_open(int n)
{
  graph_t *g;

  xnew(g, 1);
  g->n = n;
  xnew(g->mat, n * n);
  xnew(g->zhist, n + 1);
  return g;
}



static void graph_close(graph_t *g)
{
  free(g->mat);
  free(g->zhist);
  free(g);
}



static void graph_empty(graph_t *g)
{
  memset(g->mat, 0, sizeof(g->mat[0]) * g->n * g->n);
}



static void graph_link(graph_t *g, int i, int j)
{
  g->mat[i*g->n + j] = '\1';
  g->mat[j*g->n + i] = '\1';
}



__inline static int graph_linked(const graph_t *g, int i, int j)
{
  return g->mat[i*g->n + j];
}



static int graph_deg(const graph_t *g, int i)
{
  int j, deg = 0;

  for ( j = 0; j < g->n; j++ ) deg += g->mat[i*g->n + j];
  return deg;
}



static void graph_zhist(const graph_t *g)
{
  int i, n = g->n, z;

  for ( i = 0; i < n; i++ ) {
    z = graph_deg(g, i);
    g->zhist[z] += 1;
  }
}


static void graph_zhist_print(const graph_t *g)
{
  int z;

  for ( z = 0; z <= 12; z++ )
    if ( g->zhist[z] > 0 )
      printf("%2d %g\n", z, g->zhist[z]);
}



#endif /* GRAPH_H__ */

