#include "robot.h"

#include <cassert>
#include <stdexcept>

#include "context.h"
#include "../util/util.h"

using namespace arma;
using namespace std;
using namespace edu_berkeley_cs184::util;

namespace edu_berkeley_cs184 {
namespace robot {

/**
 * Takes ownership of root
 */
LinkedTreeRobot::LinkedTreeRobot(const vec3& pos, TreeNode* root) : _rootPosition(pos), _root(root) {
  assert(_root != NULL && !root->isLeafNode());

  _root->assignNodeIndicies(0, 0); // assign indicies to each of the joints and effectors separately, index based from 0.

  // compute the number of joints & number of effectors
  _numJoints    = _root->numEdges();
  _numEffectors = _root->numLeafNodes();

  _root->gatherLeaves(_effectors);
  assert(_effectors.size() == _numEffectors);
}

LinkedTreeRobot::~LinkedTreeRobot() {
  delete _root;
}

mat& LinkedTreeRobot::computeJacobian(mat& m) const {
  m.zeros(_numEffectors * 3, _numJoints);

  // for each effector, traverse up its link hierarchy and compute
  // J(i,j) = (partial S_i / partial theta_j)

  for (vector<TreeNode*>::const_iterator it = _effectors.begin();
      it != _effectors.end();
      ++it) {

    TreeNode *effector = *it;
    size_t effectorId = effector->getIdentifier();

    Context ctx(_rootPosition);
    effector->getContextForNode(ctx);

    vec3 effectorPos = ctx.getCurrentOrigin();

    TreeNode *prevNode = effector;
    TreeNode *curNode  = effector->getParent();
    while (curNode != NULL) {
      vec3 rotAxis = prevNode->getRotationAxis(ctx);
      size_t jointId = prevNode->getEdgeIdentifier();
      ctx.popContext();

      vec3 jacobianEntry = cross(rotAxis, effectorPos - ctx.getCurrentOrigin()); 
      m(3 * effectorId,     jointId) = jacobianEntry[0];
      m(3 * effectorId + 1, jointId) = jacobianEntry[1];
      m(3 * effectorId + 2, jointId) = jacobianEntry[2];

      prevNode = curNode;
      curNode  = curNode->getParent();
    }
  }

  return m;
}

vec& LinkedTreeRobot::getEffectorPositions(vec& buffer) const {
  buffer.set_size(3 * _numEffectors);
  for (size_t i = 0; i < _numEffectors; i++) {
    vec pos = _effectors[i]->getGlobalPosition(_rootPosition);
    buffer(3 * i)     = pos[0];
    buffer(3 * i + 1) = pos[1];
    buffer(3 * i + 2) = pos[2];
  }
  return buffer;
}

vector<vec3>& LinkedTreeRobot::getEffectorPositions(vector<vec3>& buf) const {
  for (vector<TreeNode*>::const_iterator it = _effectors.begin(); 
      it != _effectors.end();
      ++it) 
    buf.push_back((*it)->getGlobalPosition(_rootPosition));
  return buf;
}

static inline vec clamp(const vec& input, const double maxMag) {
  vec v = input;
  assert(v.n_elem % 3 == 0);
  for (size_t i = 0; i < v.n_elem / 3; i++) {
    vec3 cur = v.rows(3 * i, 3 * i + 2);
    if (norm(cur, 2) > maxMag) {
      normalize_vec3(cur);
      cur = maxMag * cur;
      v.rows(3 * i, 3 * i + 2) = cur;
    }
  }
  return v;
}

vec LinkedTreeRobot::computeDeltaThetas(const vec& desiredPositions) const {

  assert(_numEffectors * 3 == desiredPositions.n_elem);

  vec s(3 * _numEffectors);
  getEffectorPositions(s);

  vec e = clamp(desiredPositions - s, 1.0);

  mat J;
  computeJacobian(J);

  if (_method == PINV) { 
    //cout << "PINV" << endl;
    mat pInvJ = pinv(J); // magic. uses SVD
    vec soln = pInvJ * e;
    assert(soln.n_elem == _numJoints);
    return soln;
  } else if (_method == DLS) {
    //cout << "DLS" << endl;
    // delta theta = J^T(JJ^T + lambda^2I)^-1 * e
    
    mat JT = trans(J); 
    mat JJT = J * JT;

    mat I(J.n_rows, J.n_rows); 
    I.eye();

    double lambda = 1.5;

    double lambda_squared = lambda * lambda;

    vec soln = (JT * inv(JJT + lambda_squared * I)) * e;
    assert(soln.n_elem == _numJoints);
    return soln;

  } else
    throw runtime_error("invalid method");

}

void LinkedTreeRobot::updateThetas(const vec& deltas) {
  assert(deltas.n_elem == _numJoints);
  _root->updateThetas(deltas);
}

/** render links as lines for now */
void LinkedTreeRobot::renderRobot() const {
  Context ctx(_rootPosition);
  _root->renderTree(ctx);
}

void LinkedTreeRobot::solveIK(const vec& desired) {
  vec deltaThetas = computeDeltaThetas(desired);
  updateThetas(deltaThetas);
}

}
}
