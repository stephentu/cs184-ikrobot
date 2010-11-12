#include <stdio.h>
#include <stdlib.h>
#include <cstdarg>
#include <iostream>

#include <armadillo>

#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glut.h>

#include "robot/robot.h"
#include "robot/tree.h"
#include "util/util.h"

using namespace std;
using namespace arma;
using namespace edu_berkeley_cs184::robot;
using namespace edu_berkeley_cs184::util;

template <class T>
static inline vector<T> makeVector(size_t count, ...) {
  vector<T> vec;
  va_list ap;
  va_start(ap, count);
  for (size_t i = 0; i < count; i++) 
    vec.push_back(va_arg(ap, T));
  va_end(ap);
  return vec;
}

/** flattens a vector of vec3s into one large vec */
static inline vec flattenVector(const vector<vec3>& buffer) {
  vec v(3 * buffer.size());
  for (size_t i = 0; i < buffer.size(); i++)
    v.rows(3 * i, 3 * i + 2) = buffer[i];
  return v;
}


static LinkedTreeRobot *robot;
static vector<vec3> targets;

static int width  = 600;
static int height = 600;

static int displayMode = GL_RENDER;

#define PICK_BUFFER_SIZE 256
static unsigned int PICK_BUF[PICK_BUFFER_SIZE];

static int mouseX, mouseY; // mouseX and mouseY **AFTER** glut adjustment in y

static void setViewport() {
  glViewport(0, 0, width, height);
}

static void initScene() {
  // set background to black, fully transparent (alpha = 0)
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f); 

  // turn z-buffer depth testing on
  glEnable(GL_DEPTH_TEST);
  glClearDepth(1.0f);
  glDepthFunc(GL_LEQUAL);

  // smooth shading
  glShadeModel(GL_SMOOTH);
  
  // lighting
  GLfloat position[] = { 5.0f, 5.0f, 0.0f, 1.0f };
  glLightfv(GL_LIGHT0, GL_POSITION, position);

  GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
  GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  // material color
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  GLfloat rgba[] = { 0.8f, 0.8f, 0.8f, 1.0f };
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, rgba);
  glMateriali(GL_FRONT, GL_SHININESS, 96);

  // normal vector normalization
  glEnable(GL_NORMALIZE);
  
  // set up viewport
  setViewport();

  // pick buffer
  glSelectBuffer(PICK_BUFFER_SIZE, PICK_BUF);
  
  // set up model matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0f, 0.0f, -15.0f); // move everything into screen 15 deep
}

static void reshape(int w, int h) {
  width  = w;
  height = h;
  setViewport();
}

static void drawSphere(const vec3& center, const double radius) {
  glPushMatrix();
  glTranslated(center[0], center[1], center[2]);
  glutSolidSphere(radius, 20, 20);
  glPopMatrix();
}

/** based on: http://nehe.gamedev.net/data/articles/article.asp?article=13 */
static inline vec3 getWorldSpacePos(int x, int y, double zbuf) {
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLfloat winX, winY;
  GLdouble posX, posY, posZ;

  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);

  winX = (float)x;
  winY = (float)(viewport[3] - y);

  gluUnProject(winX, winY, zbuf, modelview, projection, viewport, &posX, &posY, &posZ);

  return makeVec3(posX, posY, posZ);
}

static void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // set up projection
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (displayMode == GL_SELECT) { // need to setup pick matrix
    GLint view[4];
    glGetIntegerv(GL_VIEWPORT, view);
    gluPickMatrix(mouseX, mouseY, 1.0, 1.0, view);
  }

  gluPerspective(45.0f, (GLfloat)width/(GLfloat)height, 1.0f, 100.0f);

  if (displayMode == GL_SELECT) { // init name stack
    glInitNames();
    glPushName(0); // need at least 1 element on stack
  }

  glMatrixMode(GL_MODELVIEW);
  if (displayMode == GL_SELECT) { // only render set points
    for (size_t i = 0; i < targets.size(); i++) {
      glLoadName(i);
      glColor3f(1.0, 0.0, 0.0);
      drawSphere(targets[i], 0.1);
    }
  } else { // render everything
    for (size_t i = 0; i < targets.size(); i++) {
      glColor3f(1.0, 0.0, 0.0);
      drawSphere(targets[i], 0.1);
    }
    robot->renderRobot();
  }

  if (displayMode == GL_RENDER)
    glutSwapBuffers();
}

static void redisplay() {
  glutPostRedisplay();
}

static void specialKeyboardHandler(int key, int x, int y) {
  int mod = glutGetModifiers();
  if (mod == GLUT_ACTIVE_SHIFT) {
    switch (key) {
      case GLUT_KEY_LEFT: // -x translate
        glMatrixMode(GL_MODELVIEW);
        glTranslated(-0.5, 0.0, 0.0);
        break;
      case GLUT_KEY_RIGHT: // +x translate
        glMatrixMode(GL_MODELVIEW);
        glTranslated(0.5, 0.0, 0.0);
        break;
      case GLUT_KEY_UP: // -y translate
        glMatrixMode(GL_MODELVIEW);
        glTranslated(0.0, -0.5, 0.0);
        break;
      case GLUT_KEY_DOWN: // +y translate
        glMatrixMode(GL_MODELVIEW);
        glTranslated(0.0, 0.5, 0.0);
        break;
    }
  } else {
    switch (key) {
      case GLUT_KEY_LEFT: // -y axis rot
        glMatrixMode(GL_MODELVIEW);
        glRotated(-10.0, 0.0, 1.0, 0.0);
        break;
      case GLUT_KEY_RIGHT: // +y axis rot
        glMatrixMode(GL_MODELVIEW);
        glRotated(10.0, 0.0, 1.0, 0.0);
        break;
      case GLUT_KEY_UP: // -x axis rot
        glMatrixMode(GL_MODELVIEW);
        glRotated(-10.0, 1.0, 0.0, 0.0);
        break;
      case GLUT_KEY_DOWN: // +x axis rot
        glMatrixMode(GL_MODELVIEW);
        glRotated(10.0, 1.0, 0.0, 0.0);
        break;
    }
  }
}

static void keyboardHandler(unsigned char key, int x, int y) {
  switch (key) {
    case '+': // zoom in
      glMatrixMode(GL_MODELVIEW);
      glScaled(1.5, 1.5, 1.5);
      break;
    case '-': // zoom out
      glMatrixMode(GL_MODELVIEW);
      double two_thirds = 2.0 / 3.0;
      glScaled(two_thirds, two_thirds, two_thirds);
      break;
  }
}

static int lastMouseButton, lastMouseState;

#define LEFT_DOWN_CLICK (lastMouseButton == GLUT_LEFT_BUTTON && lastMouseState == GLUT_DOWN)
#define LEFT_UP_CLICK (lastMouseButton == GLUT_LEFT_BUTTON && lastMouseState == GLUT_UP)

static int activeTarget = -1;
static float activeTargetZBuf = 0.0;

static void mouseClicked(int button, int state, int x, int y) {
  lastMouseButton = button;
  lastMouseState  = state;

  mouseX = x;
  mouseY = height - y;

  if (LEFT_DOWN_CLICK) {
    displayMode = GL_SELECT;
    glRenderMode(GL_SELECT);
    display();
    displayMode = GL_RENDER;
    int numHits = glRenderMode(GL_RENDER);
    cout << "numHits: " << numHits << endl;

    if (numHits > 0) { // find the item which has the smallest zmin
      assert(PICK_BUF[0] == 1);
      unsigned int zminMinSoFar = PICK_BUF[1];
      unsigned int minIdx = 0;

      for (int i = 1; i < numHits; i++) {
        assert(PICK_BUF[4 * i] == 1);
        if (PICK_BUF[4 * i + 1] < zminMinSoFar) {
          zminMinSoFar = PICK_BUF[4 * i + 1]; 
          minIdx = i;
        }
      }

      unsigned int item = PICK_BUF[4 * minIdx + 3];
      activeTarget = item;
      glReadPixels(mouseX, mouseY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &activeTargetZBuf);
    }
    redisplay();
  } else if (LEFT_UP_CLICK) { // reset the active cube
    activeTarget = -1;
  }
}

static void mouseDragged(int x, int y) {
  if (LEFT_DOWN_CLICK && activeTarget != -1) {
    //cout << "moueDragged: " << x << ", " << y << endl;
    vec3 pos = getWorldSpacePos(x, y, activeTargetZBuf);
    //cout << pos << endl;
    targets[activeTarget] = pos;

    robot->solveIK(flattenVector(targets));
  }
}


int main(int argc, char **argv) {  

  TreeNode* root = new INode(
    makeVector<LinkState*>(2, new LinkState(1.0, makeVec3(0, 0, 1), 0.0, makeVec3(-1, -1, 0)), new LinkState(1.0, makeVec3(0, 0, 1), 0.0, makeVec3(1, -1, 0))), 
    makeVector<TreeNode*>(2, new LNode(), new LNode()));

  robot = new LinkedTreeRobot(makeVec3(0, 0, 0), root);

  cout << "numJoints: " << robot->getNumJoints() << ", " <<
          "numEffectors: " << robot->getNumEffectors() << endl;

  robot->getEffectorPositions(targets);

  // -----------------

  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(width, height);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("CS184 Final Project - Inverse Kinematics");

  initScene();

  glutIdleFunc(redisplay);
  glutKeyboardFunc(keyboardHandler);
  glutSpecialFunc(specialKeyboardHandler);
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutMouseFunc(mouseClicked);
  glutMotionFunc(mouseDragged);

  glutMainLoop();

  return 0;
}
