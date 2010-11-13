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

pair<size_t, size_t> INode::assignNodeIndicies(const size_t iNodeOffset, const size_t lNodeOffset) {
  // assign ids to all the children first
  size_t curOffset = 0;
  for (size_t i = 0; i < _states.size(); i++) 
    curOffset += _states[i]->assignJointIndicies(iNodeOffset + curOffset);
  
  // now do this recursively with bookkeeping
  size_t numINodes = curOffset;
  size_t numLNodes = 0;
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it) {
    pair<size_t, size_t> thisNode = (*it)->assignNodeIndicies(iNodeOffset + numINodes, 
                                                              lNodeOffset + numLNodes);
    numINodes += thisNode.first;
    numLNodes += thisNode.second;
  }
  return pair<size_t, size_t>(numINodes, numLNodes);
}

vector<TreeNode*>& INode::gatherLeaves(vector<TreeNode*>& buffer) {
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    (*it)->gatherLeaves(buffer);
  return buffer;
}

size_t INode::getIdentifier() const {
  throw runtime_error("getIdentifier on intermediate node");
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
    //glBegin(GL_LINES);
    //glColor3d(0.0, 1.0, 0.0);
    //glVertex3d(pos[0], pos[1], pos[2]);
    //glVertex3d(endpoint[0], endpoint[1], endpoint[2]);
    //glEnd();

    // calculate orthonormal basis for cylinder on joint
    vec3 u, v, n;
    _states[i]->getBasis(ctx, u, v, n);

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

    //cout << "--" << endl;
    //for (int iii = 0; iii < 16; iii++) {
    //  cout << m[iii] << endl;
    //}
    //cout << "--" << endl;

    if (isRootNode()) {
      glColor3d(0.0, 0.0, 0.8);
      glPushMatrix();
        glTranslated(startpoint[0], startpoint[1], startpoint[2]);
        glutSolidSphere(0.1, 20, 20);
      glPopMatrix();
    }

    GLUquadricObj *quadric = gluNewQuadric();

    glPushMatrix();
      glColor3d(0.0, 1.0, 0.0);
      glTranslated(startpoint[0], startpoint[1], startpoint[2]);
      glMultMatrixd(m);
      gluCylinder(quadric, 0.05, 0.05, _states[i]->getLength(), 20, 20);
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

pair<size_t, size_t> LNode::assignNodeIndicies(const size_t iNodeOffset, const size_t lNodeOffset) {
  id = lNodeOffset;
  return pair<size_t, size_t>(0, 1);
}

vector<TreeNode*>& LNode::gatherLeaves(vector<TreeNode*>& buffer) {
  buffer.push_back(this);
  return buffer;
}

size_t LNode::getIdentifier() const { return id; }

void LNode::updateThetas(const vec& deltas, const vec& axes) {}
void LNode::renderTree(Context& ctx) const {}

}
}
