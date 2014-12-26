/* shape sampler */



#include "mtrand.h" /* Mersenne-Twister random number generator */
#include "util.h" /* utility functions */
#include "shaper.h"



/* parameters for the mesh */

/* theta */
#define NPART 45
#define NSUB  4
#define N     (NPART * NSUB)
/* phi */
#define MPART 15
#define MSUB  6
#define M     (MPART * MSUB)

#define PX    4
#define PY    4

int stop = 0;
int interval = 50;
int printtext = 1;

int lmax = 7; /* maximal l of spherical harmonics */



GLfloat points[N+1][M+1][3];
shaper_t *shaper;



/* calculate the coordinates of the points */
static void genpts(GLfloat *x, int n, int m)
{
  int i, j, id;
  double theta, sth, cth, phi, r;

  id = 0;
  for (i = 0; i <= n; i++) { /* theta */
    theta = PI*i/n;
    sth = sin(theta);
    cth = cos(theta);
    shaper_getal(shaper, cth);

    for (j = 0; j <= m; j++) { /* phi */
      phi = 2*PI*j/m;
      r = shaper_getr(shaper, phi, 1, 0.3);
      x[id*3+0] = (GLfloat)( r * sth * cos(phi) );
      x[id*3+1] = (GLfloat)( r * sth * sin(phi) );
      x[id*3+2] = (GLfloat)( r * cth );
      id++;
    }
  }
}



void timer(int value)
{
  glutPostRedisplay(); /* tell window to redraw */
  glutTimerFunc(interval, timer, ++value);
}



void display(void)
{
  static double t = 0, dt = 0.01;
  int i, j;
  double scl;

  if (stop) return;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  t += dt;

  /* sample */
  //for ( i = 0; i < 10; i++ ) shaper_mc(shaper, 0.5);
  shaper_lang(shaper, 0.1);

  /* generate points */
  genpts((GLfloat *) points, N, M);

  glPushMatrix();
  glRotatef(-60.0, 1, 0, 0);
  glRotatef(-t*300, sin(t*2), 0, cos(t*2)); /* view point */
  scl = 0.8;
  glScaled(scl, scl, scl);

  for (i = 0; i < NPART; i++) {
    GLfloat color[4] = {1.f, 1.5f, 0, 0};
    huetorgb(color, 1.f*i/NPART);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, color);
    glLightfv(GL_LIGHT0, GL_AMBIENT, color);
    for (j = 0; j < MPART; j++) {
      glMap2f(GL_MAP2_VERTEX_3,
          0.0f, 1.0f, /* u parameter for phi */
          3, /* three floats per phi value in the data array */
          MSUB + 1, /* number of control points (order) along u */
          0.0f, 1.0f, /* v parameter for theta */
          3*(M + 1), /* 3*(M+1) floats per theta value */
          NSUB + 1,  /* number of control points (order) along v */
          &points[NSUB*i][MSUB*j][0]); /* address of this block */
      /* create a 2D grid of PX x PY points in [0.0, 1.0] x [0.0, 1.0] */
      glMapGrid2f(PX, 0.0, 1.0, PY, 0.0, 1.0);
      /* evaluate the grid points */
      glEvalMesh2(GL_FILL, 0, PX, 0, PY);
    }
  }
  glPopMatrix();

  ///* print text */
  //if ( printtext ) {
  //  char s[80], *p;
  //  sprintf(s, "Y(%2d, %2d) -> Y(%2d, %2d), %3.f%%",
  //      l0, m0, l1, m1, 100.0*t/(2*PI));
  //  glDisable(GL_LIGHTING);
  //  glColor3f(0.5, 0.5, 0.5);
  //  glRasterPos2f(-0.98, -0.9);
  //  for (p = s; *p != '\0'; p++) glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
  //  glEnable(GL_LIGHTING);
  //}

  glutSwapBuffers();
}



static void init(void)
{
  GLfloat position[] = {-.6f, 0, -1, 0};

  shaper = shaper_open(lmax);

  srand(time(NULL));
  glClearColor(0.f, 0.f, 0.f, 0.f); /* use black to clear screen */

  glLightfv(GL_LIGHT0, GL_POSITION, position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MAP2_VERTEX_3);
  glEnable(GL_AUTO_NORMAL);

  glutTimerFunc(interval, timer, 0);
}



static void done(void)
{
  shaper_close(shaper);
}



static void reshape(int w, int h)
{
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.0f*w/h, 1.0f*w/h, -1.0f, 1.0f, -5.0f, 5.0f);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}



static void keyboard(unsigned char key, int x, int y)
{
  (void) x; (void) y;
  if (key == 27) {
    exit(0);
  } if (key == ' ') {
    stop = !stop;
    if (!stop) timer(0);
  } else if (key == 'f') {
    static int full = 0, x0, y0, w, h;
    full = !full;
    if (full) {
      x0 = glutGet(GLUT_WINDOW_X);
      y0 = glutGet(GLUT_WINDOW_Y);
      w = glutGet(GLUT_WINDOW_WIDTH);
      h = glutGet(GLUT_WINDOW_HEIGHT);
      glutFullScreen();
    } else {
      glutPositionWindow(x0, y0);
      glutReshapeWindow(w, h);
    }
  }
}



int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);
  init();
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMainLoop();
  done();
  return 0;
}

