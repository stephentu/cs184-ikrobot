#include <stdio.h>
#include <stdlib.h>

#include <cstdarg>
#include <iostream>
#include <cmath>
#include <ctime>

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

/** http://cc.byexamples.com/2007/05/25/nanosleep-is-better-than-sleep-and-usleep */
static int msleep(unsigned long milisec) {
    struct timespec req={0};
    time_t sec=(int)(milisec/1000);
    milisec=milisec-(sec*1000);
    req.tv_sec=sec;
    req.tv_nsec=milisec*1000000L;
    while(nanosleep(&req,&req)==-1)
         continue;
    return 1;
}

static LinkedTreeRobot *robot;
static vector<vec3> targets;

static int width  = 600;
static int height = 600;

static int displayMode = GL_RENDER;

#define PICK_BUFFER_SIZE 256
static unsigned int PICK_BUF[PICK_BUFFER_SIZE];

static int mouseX, mouseY; // mouseX and mouseY **AFTER** glut adjustment in y

static int solnType = 0;

static bool isPlayback = false; // if false, in rec mode. if true, in playback mode

static size_t playbackFps = 20; // 20 frames per second

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

    vector<vec3> nodePositions;
    robot->getNodePositions(nodePositions);

    for (size_t i = 0; i < nodePositions.size(); i++) {
      glLoadName(targets.size() + i);
      glColor3f(1.0, 0.0, 0.0);
      drawSphere(nodePositions[i], 0.08);
    }
  } else { // render everything
    if (!isPlayback)
      for (size_t i = 0; i < targets.size(); i++) {
        glColor3f(1.0, 0.0, 0.0);
        drawSphere(targets[i], 0.1);
      }
    else if (isPlayback && !robot->isEndOfBuffer()) { // if we are in playback and there are contents in the buffer
      robot->restore(); // restore robot to next animation frame
      robot->getEffectorPositions(targets); // reset effectors
      robot->advance(); // push frame pointer forward
      if (robot->isEndOfBuffer()) robot->resetPointer(); // wrap around if necessary
    }
    robot->renderRobot(); // draw the robot to the buffer
    if (isPlayback) { // impose animation time delay
      const double secPerFrame = 1.0 / ((double)playbackFps);
      msleep((long)(secPerFrame * 1000.0));
    }
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

static const double two_thirds = 2.0 / 3.0;

static void keyboardHandler(unsigned char key, int x, int y) {
  switch (key) {
    case '+': // zoom in
      glMatrixMode(GL_MODELVIEW);
      glScaled(1.5, 1.5, 1.5);
      break;
    case '-': // zoom out
      glMatrixMode(GL_MODELVIEW);

      glScaled(two_thirds, two_thirds, two_thirds);
      break;
    case 's': // switch soln type
      solnType = (solnType + 1) % ((int)NUM_TYPES);
      robot->setMethod((SolnType)solnType);
      break;
    case 'm': // mark
      robot->mark();
      break;
    case 'r': // reset
      robot->reset();
      robot->getEffectorPositions(targets); // also reset effectors
      break;
    case 'a': // append
      robot->save();
      break;
    case 'x': // next
      isPlayback = !isPlayback;
      break;
  }
}

static int lastMouseButton, lastMouseState;

#define LEFT_DOWN_CLICK (lastMouseButton == GLUT_LEFT_BUTTON && lastMouseState == GLUT_DOWN)
#define LEFT_UP_CLICK (lastMouseButton == GLUT_LEFT_BUTTON && lastMouseState == GLUT_UP)

static bool hasActiveTarget = false;
static unsigned int activeTarget = 0;
static float activeTargetZBuf = 0.0;

static void innerNodeClicked(size_t idx, int x, int y) {
  cout << "inner node: " << (activeTarget-targets.size()) << endl;
  robot->toggleConstraint(idx);
}

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
      hasActiveTarget = true;
      glReadPixels(mouseX, mouseY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &activeTargetZBuf);

      if (activeTarget >= targets.size())
        innerNodeClicked(activeTarget - targets.size(), x, y);
    }
    redisplay();
  } else if (LEFT_UP_CLICK) { // reset the active cube
    hasActiveTarget = false;
  }
}

static void mouseDragged(int x, int y) {
  if (LEFT_DOWN_CLICK && hasActiveTarget) {

    if (activeTarget < targets.size()) {
      //cout << "moueDragged: " << x << ", " << y << endl;
      //cout << pos << endl;
      vec3 pos = getWorldSpacePos(x, y, activeTargetZBuf);
      targets[activeTarget] = pos;
      robot->solveIKWithConstraints(flattenVector(targets));
      //robot->solveIK(flattenVector(targets));
    } else if (activeTarget == targets.size()) { // HACK: the root inner node is always the first node after the effectors
      vec3 pos = getWorldSpacePos(x, y, activeTargetZBuf);
      robot->setRootPosition(pos);
      robot->solveIKWithConstraints(flattenVector(targets));
    }
  }
}


int main(int argc, char **argv) {  

  // 3 link model
  //TreeNode* root = new INode(
  //  makeVector<LinkState*>(1, new LinkState(1.0, makeVec3(0, 0, 1), makeVec3(-1, 0, 0))),
  //  makeVector<TreeNode*>(1, 
  //    new INode(
  //      makeVector<LinkState*>(1, new LinkState(2.0, makeVec3(-1, 0, 0), makeVec3(0, -1, 0))),
  //      makeVector<TreeNode*>(1, 
  //        new INode(
  //          makeVector<LinkState*>(1, new LinkState(1.0, makeVec3(0, -1, 0), makeVec3(0, 0, 1))),
  //          makeVector<TreeNode*>(1, new LNode()))))));

  vec3 zhat = makeVec3(0, 0, 1);
  
  // spider
  //TreeNode *root = new INode(
  //  makeVector<LinkState*>(2, new LinkState(1, zhat, makeVec3(sqrt(3.0)/2.0, -0.5, 0)), new LinkState(1, zhat, makeVec3(-sqrt(3.0)/2.0, -0.5, 0))),
  //  makeVector<TreeNode*>(2, 
  //    new INode(
  //      makeVector<LinkState*>(1, new LinkState(1, zhat, makeVec3(0.5, -sqrt(3.0)/2.0, 0))),
  //      makeVector<TreeNode*>(1, new LNode())),
  //    new INode(
  //      makeVector<LinkState*>(1, new LinkState(1, zhat, makeVec3(-0.5, -sqrt(3.0)/2.0, 0))),
  //      makeVector<TreeNode*>(1, new LNode()))));

  // ball joint
  //TreeNode *root = new INode(
  //  makeVector<LinkState*>(1, new AxisBallAndSocketJoint(1.0, makeVec3(0, 1, 0))),
  //  makeVector<TreeNode*>(1, 
  //    //new INode(
  //    //  makeVector<LinkState*>(1, new RotationJoint(1.0, zhat, makeVec3(0, 1, 0))),
  //    //  makeVector<TreeNode*>(1, new LNode())
  //    //)
  //    new LNode()
  //  )
  //);

  //TreeNode *root = new INode(
  //  makeVector<LinkState*>(1, new RotationJoint(1.0, zhat, makeVec3(1, 0, 0))),
  //  makeVector<TreeNode*>(1, 
  //    new INode( 
  //      makeVector<LinkState*>(1, new RotationJoint(1.0, zhat, makeVec3(1, 0, 0))),
  //      makeVector<TreeNode*>(1, new LNode())
  //    
  //    )));

  // humanoid

  TreeNode *root = new INode(
    makeVector<LinkState*>(4, 
      new AxisBallAndSocketJoint(1.0, makeVec3(1, 0, 0)),
      new AxisBallAndSocketJoint(0.5, makeVec3(0, 1, 0)),
      new AxisBallAndSocketJoint(1.0, makeVec3(-1, 0, 0)),
      new AxisBallAndSocketJoint(2.0, makeVec3(0, -1, 0))
    ),
    makeVector<TreeNode*>(4,
      new INode(
        makeVector<LinkState*>(1, new RotationJoint(0.8, zhat, makeVec3(0, -1, 0))),
        makeVector<TreeNode*>(1, new LNode())
      ),
      new LNode(),
      new INode(
        makeVector<LinkState*>(1, new RotationJoint(0.8, zhat, makeVec3(0, -1, 0))),
        makeVector<TreeNode*>(1, new LNode())
      ),
      new INode(
        makeVector<LinkState*>(2, 
          new AxisBallAndSocketJoint(1.2, makeVec3(1, 0, 0)),
          new AxisBallAndSocketJoint(1.2, makeVec3(-1, 0, 0))
        ),
        makeVector<TreeNode*>(2, 
          new INode(
            makeVector<LinkState*>(1, new RotationJoint(0.8, zhat, makeVec3(0, -1, 0))),
            makeVector<TreeNode*>(1, new LNode())
          ),
          new INode(
            makeVector<LinkState*>(1, new RotationJoint(0.8, zhat, makeVec3(0, -1, 0))),
            makeVector<TreeNode*>(1, new LNode())
          )
        )
      )
    )
  );

  robot = new LinkedTreeRobot(makeVec3(0, 0, 0), root);

  cout << "numJoints: " << robot->getNumJoints() << ", " <<
          "numEffectors: " << robot->getNumEffectors() << endl;

  robot->getEffectorPositions(targets);

  cout << "Got effector positions" << endl;

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
