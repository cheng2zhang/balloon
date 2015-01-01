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


int lmax = 9; /* maximal l of spherical harmonics */
double dt = 0.01;
double tp = 1;

enum {
  MENU_FIRST = 1000,
  MENU_SCALE_DOWN,
  MENU_SCALE_UP,
  MENU_FULLSCREEN,
  MENU_LAST
};

int mouse_x, mouse_y, mouse_down;
int stop = 0;
int interval = 50;
int printtext = 1;



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
  int i, j;
  static int once;

  if ( stop ) return;

  if ( !once ) {
    GLfloat s = 0.7f;
    glRotatef(-60.f, 1.f, 0.f, 0.f);
    glScalef(s, s, s);
    once = 1;
  }

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* sample */
  //for ( i = 0; i < 10; i++ ) shaper_mc(shaper, 0.5);
  shaper_lang(shaper, dt, tp);

  /* generate points */
  genpts((GLfloat *) points, N, M);

  for (i = 0; i < NPART; i++) {
    GLfloat color[4];
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

  glutSwapBuffers();
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



/* toggle the full screen state */
static void fullscreen(void)
{
  static int full = 0, x, y, w, h;

  full = !full;
  if (full) {
    x = glutGet(GLUT_WINDOW_X);
    y = glutGet(GLUT_WINDOW_Y);
    w = glutGet(GLUT_WINDOW_WIDTH);
    h = glutGet(GLUT_WINDOW_HEIGHT);
    glutFullScreen();
  } else {
    glutPositionWindow(x, y);
    glutReshapeWindow(w, h);
  }
}



/* menu function */
static void menu(int id)
{
  if ( id == MENU_SCALE_DOWN || id == MENU_SCALE_UP ) { /* scaling */
    GLfloat s = (id == MENU_SCALE_UP) ? 1.05f : 1.0f/1.05f;
    glScalef(s, s, s);
    glutPostRedisplay();
  } else if ( id == MENU_FULLSCREEN ) { /* full screen */
    fullscreen();
  }
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
    fullscreen();
  }
}



static void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    mouse_down++;
    if (button == 3)
      menu(MENU_SCALE_UP);
    else if (button == 4)
      menu(MENU_SCALE_DOWN);
  } else if (--mouse_down <= 0) {
    mouse_down = 0;
  }
  mouse_x = x;
  mouse_y = y;
}



/* mouse motion function for GLUT */
static void motion(int x, int y)
{
  if ( x == mouse_x && y == mouse_y )
    return;

  if ( mouse_down ) {
    float angx = (float)( (y - mouse_y) * 360.f / glutGet(GLUT_WINDOW_HEIGHT) );
    float angy = (float)( (x - mouse_x) * 360.f / glutGet(GLUT_WINDOW_WIDTH) );
    float mat[4][4];

    glGetFloatv(GL_MODELVIEW_MATRIX, (GLfloat *) mat);
    glRotatef(angx, mat[0][0], mat[1][0], mat[2][0]);
    glRotatef(angy, mat[0][1], mat[1][1], mat[2][1]);
    glutPostRedisplay();
  }
  mouse_x = x;
  mouse_y = y;
}



static void init(void)
{
  GLfloat position[] = {-.6f, 0, -1, 0};

  shaper = shaper_open(lmax);

  glClearColor(0.f, 0.f, 0.f, 0.f); /* use black to clear screen */

  glLightfv(GL_LIGHT0, GL_POSITION, position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MAP2_VERTEX_3);
  glEnable(GL_AUTO_NORMAL);

  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutCreateMenu(menu);
  glutAddMenuEntry("Shrink",      MENU_SCALE_DOWN);
  glutAddMenuEntry("Enlarge",     MENU_SCALE_UP);
  glutAddMenuEntry("FullScreen",  MENU_FULLSCREEN);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  glutTimerFunc(interval, timer, 0);
}



static void done(void)
{
  shaper_close(shaper);
}



int main(int argc, char **argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(600, 600);
  glutInitWindowPosition(100, 100);
  glutCreateWindow(argv[0]);
  init();
  glutMainLoop();
  done();
  return 0;
}

