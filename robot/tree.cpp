#include "tree.h"

#include <stdexcept>
#include <cassert>
#include <stack>

#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glut.h>

#include "../util/util.h"

using namespace std;
using namespace arma;
using namespace edu_berkeley_cs184::util;

namespace edu_berkeley_cs184 {
namespace robot {

TreeNode::TreeNode() { _parent = NULL; }
TreeNode::~TreeNode() {}

vec3 TreeNode::getGlobalPosition(const vec3& rootPoint) {
  Context ctx(rootPoint);
  getContextForNode(ctx);
  return ctx.getCurrentOrigin();
}

Context& TreeNode::getContextForNode(Context& ctx) {
  stack<TreeNode*> pathToRoot;
  stack<size_t> indexes;
  TreeNode* cur = getParent();
  TreeNode* prev = this;
  while (cur != NULL) {
    pathToRoot.push(cur);
    indexes.push(prev->getIndex());
    prev = cur;
    cur = cur->getParent();
  }

  assert(pathToRoot.size() == indexes.size());

  while (!pathToRoot.empty()) {
    TreeNode *n = pathToRoot.top();
    size_t idx = indexes.top();
    pathToRoot.pop(); indexes.pop();
    vector<LinkState*>::const_iterator state_iter = n->getLinkStates() + idx;
    LinkState* state = *state_iter;
    state->pushContext(ctx);
  }

  return ctx;
}


/**
 * Takes ownership of objs in states and objs in children
 */
INode::INode(const vector<LinkState*>& states, 
             const vector<TreeNode*>& children) : _states(states), _kids(children) {
  assert(states.size() == children.size());
  assert(states.size() > 0); // INodes are required to have >= 1 child
  for (size_t i = 0; i < _kids.size(); i++) {
    _kids[i]->setParent(this);
    _kids[i]->setIndex(i);
  }
}

/**
 * Recursively deletes link states first, and then children nodes
 */
INode::~INode() {
  for (vector<LinkState*>::iterator it = _states.begin();
      it != _states.end();
      ++it)
    delete *it;
  for (vector<TreeNode*>::iterator it = _kids.begin();
      it != _kids.end();
      ++it)
    delete *it;
}

bool INode::isLeafNode() const { return false; }

vector<TreeNode*>::const_iterator INode::getChildren() const { return _kids.begin(); }
vector<LinkState*>::const_iterator INode::getLinkStates() const { return _states.begin(); }

size_t INode::numLeafNodes() const {
  size_t sumSoFar = 0;
  for (vector<TreeNode*>::const_iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    sumSoFar += (*it)->numLeafNodes();
  return sumSoFar;
}

size_t INode::numEdges() const {
  size_t sumSoFar = _kids.size();
  for (vector<TreeNode*>::const_iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    sumSoFar += (*it)->numEdges();
  return sumSoFar;
}

size_t INode::numDOF() const {
  size_t sumSoFar = 0;
  for (vector<LinkState*>::const_iterator it = _states.begin();
      it != _states.end();
      ++it)
    sumSoFar += (*it)->dof();
  for (vector<TreeNode*>::const_iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    sumSoFar += (*it)->numDOF();
  return sumSoFar;
}

void INode::assignNodeIndicies(const size_t innerNodeOffset, 
                               const size_t jointNodeOffset,
                               const size_t leafNodeOffset,
                               size_t& numInnerNodes,
                               size_t& numJointNodes,
                               size_t& numLeafNodes) {
  numLeafNodes = 0;

  // take an id
  numInnerNodes = 1;
  identity = innerNodeOffset;

  // assign ids to all the joints first
  numJointNodes = 0;
  for (size_t i = 0; i < _states.size(); i++) 
    numJointNodes += _states[i]->assignJointIndicies(jointNodeOffset + numJointNodes);
  
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it) {

    size_t childNumInnerNodes;
    size_t childNumJointNodes;
    size_t childNumLeafNodes;

    (*it)->assignNodeIndicies(innerNodeOffset + numInnerNodes, 
                              jointNodeOffset + numJointNodes, 
                              leafNodeOffset + numLeafNodes, 
                              childNumInnerNodes,
                              childNumJointNodes,
                              childNumLeafNodes);

    numInnerNodes += childNumInnerNodes;
    numJointNodes += childNumJointNodes;
    numLeafNodes  += childNumLeafNodes;
  }
}

vector<TreeNode*>& INode::gatherLeaves(vector<TreeNode*>& buffer) {
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    (*it)->gatherLeaves(buffer);
  return buffer;
}

vector<TreeNode*>& INode::gatherInnerNodes(vector<TreeNode*>& buffer) {
  buffer.push_back(this);
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    (*it)->gatherInnerNodes(buffer);
  return buffer;
}

void INode::updateThetas(const vec& deltas, const vec& axes) {
  for (size_t i = 0; i < _states.size(); i++)
    _states[i]->updateThetas(deltas, axes);
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    (*it)->updateThetas(deltas, axes);
}

void INode::renderTree(Context& ctx) const {
  for (size_t i = 0; i < _states.size(); i++) {
    // start point
    vec3 startpoint = ctx.getCurrentOrigin();

    // calculate end point of joint
    vec3 endpoint = _states[i]->getEndpoint(ctx);

    // draw the link
    glBegin(GL_LINES);
      glColor3d(1.0, 1.0, 1.0);
      glVertex3d(startpoint[0], startpoint[1], startpoint[2]);
      glVertex3d(endpoint[0], endpoint[1], endpoint[2]);
    glEnd();

    // calculate orthonormal basis for cylinder on joint
    vec3 u, v, n;
    _states[i]->getBasis(ctx, u, v, n);

    // check if basis is really orthonormal
    assert(double_equals(dot(u, v), 0));
    assert(double_equals(dot(u, n), 0));
    assert(double_equals(dot(v, n), 0));

    assert(double_equals(norm(u, 2), 1));
    assert(double_equals(norm(v, 2), 1));
    assert(double_equals(norm(n, 2), 1));

    //cout << "pos:" << endl << pos << endl;

    //cout << "u:" << endl << u << endl;
    //cout << "v:" << endl << v << endl;
    //cout << "n:" << endl << n << endl;

    vec3 x = makeVec3(1, 0, 0); 
    vec3 y = makeVec3(0, 1, 0);
    vec3 z = makeVec3(0, 0, 1);

    double ux = dot(x, u);
    double uy = dot(y, u);
    double uz = dot(z, u);

    double vx = dot(x, v);
    double vy = dot(y, v);
    double vz = dot(z, v);

    double nx = dot(x, n);
    double ny = dot(y, n);
    double nz = dot(z, n);

    // change of orthonormal basis from uvn -> xyz
    GLdouble m[16];
    m[0]  = ux;
    m[1]  = uy;
    m[2]  = uz;
    m[3]  = 0;
    
    m[4]  = vx; 
    m[5]  = vy;
    m[6]  = vz;
    m[7]  = 0;

    m[8]  = nx;
    m[9]  = ny;
    m[10] = nz;
    m[11] = 0;

    m[12] = 0; m[13] = 0; m[14] = 0; m[15] = 1;

    mat44 A; 
    A << ux << vx << nx << 0 << endr 
      << uy << vy << ny << 0 << endr
      << uz << vz << nz << 0 << endr 
      << 0  << 0  << 0  << 1 << endr;

    //if (!double_equals(det(A), 1))
    //  cout << "A is: " << endl << A << endl;

    //cout << "det(A): " << det(A) << endl;
    const double dA = det(A);
    if (!double_equals(dA, 1)) {
      cerr << "ERROR: det(A) = " << dA << endl; 
      throw runtime_error("determinant not 1 for rotation matrix");
    }


    //cout << "--" << endl;
    //for (int iii = 0; iii < 16; iii++) {
    //  cout << m[iii] << endl;
    //}
    //cout << "--" << endl;

    if (isRootNode())
      glColor3d(0.0, 0.0, 0.8);
    else
      glColor3d(0.0, 0.0, 1.0);

    glPushMatrix();
      glTranslated(startpoint[0], startpoint[1], startpoint[2]);
      if (isRootNode())
        glutSolidSphere(0.1, 20, 20);
      else
        glutSolidSphere(0.08, 20, 20);
    glPopMatrix();

    GLUquadricObj *quadric = gluNewQuadric();

    glPushMatrix();
      glColor3d(0.0, 1.0, 0.0);
      glTranslated(startpoint[0], startpoint[1], startpoint[2]);
      glMultMatrixd(m);
      gluCylinder(quadric, 0.05, 0.05, _states[i]->getLength(), 32, 32);
    glPopMatrix();

    gluDeleteQuadric(quadric);

    // recurse into child
    _states[i]->pushContext(ctx);
      _kids[i]->renderTree(ctx);
    ctx.popContext();
  }
}

bool LNode::isLeafNode() const { return true; }

vector<TreeNode*>::const_iterator LNode::getChildren() const { 
  throw runtime_error("getChildren on leaf node"); 
}

vector<LinkState*>::const_iterator LNode::getLinkStates() const { 
  throw runtime_error("getLinkStates on leaf node"); 
}

size_t LNode::numLeafNodes() const { return 1; }

size_t LNode::numEdges() const { return 0; }

size_t LNode::numDOF() const { return 0; }

void LNode::assignNodeIndicies(const size_t innerNodeOffset, 
                               const size_t jointOffset,
                               const size_t leafNodeOffset,
                               size_t& numInnerNodes,
                               size_t& numJointNodes,
                               size_t& numLeafNodes) {
  identity = leafNodeOffset;
  numInnerNodes = 0;
  numJointNodes = 0;
  numLeafNodes  = 1;
}

vector<TreeNode*>& LNode::gatherLeaves(vector<TreeNode*>& buffer) {
  buffer.push_back(this);
  return buffer;
}

vector<TreeNode*>& LNode::gatherInnerNodes(vector<TreeNode*>& buffer) {
  return buffer;
}

void LNode::updateThetas(const vec& deltas, const vec& axes) {}
void LNode::renderTree(Context& ctx) const {}

}
}
