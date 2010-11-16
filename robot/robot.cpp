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

  size_t nInnerNodes, nJoints, nLeaves;
  _root->assignNodeIndicies(0, 0, 0, nInnerNodes, nJoints, nLeaves); 

  cout << "nInnerNodes: " << nInnerNodes << endl
       << "nJoints: " << nJoints << endl
       << "nLeaves: " << nLeaves << endl;

  // compute the number of joints & number of effectors
  _numJoints    = _root->numDOF();
  _numEffectors = _root->numLeafNodes();

  _root->gatherLeaves(_effectors);
  assert(_effectors.size() == _numEffectors);

  _root->gatherInnerNodes(_innerNodes);

  assert(nJoints == _numJoints);
  assert(nLeaves == _numEffectors);
  assert(nInnerNodes == _innerNodes.size());
}

LinkedTreeRobot::~LinkedTreeRobot() {
  delete _root;
}

mat& LinkedTreeRobot::computeJacobian(const vec& desired, mat& m, vec& axes) const {
  //cout << "--- computing jacobian ---" << endl;

  assert(desired.n_elem == 3 * _numEffectors);
  m.zeros(_numEffectors * 3, _numJoints);
  axes.zeros(_numJoints * 3);

  // for each effector, traverse up its link hierarchy and compute
  // J(i,j) = (partial S_i / partial theta_j)

  for (vector<TreeNode*>::const_iterator it = _effectors.begin();
      it != _effectors.end();
      ++it) {

    TreeNode *effector = *it;
    size_t effectorId = effector->getIdentifier();

    //cout << "considering effector: " << effectorId << endl;

    Context ctx(_rootPosition);
    effector->getContextForNode(ctx);

    vec3 effectorPos = ctx.getCurrentOrigin();

    //cout << "effector currently located at: " << effectorPos << endl;

    TreeNode *prevNode = effector;
    TreeNode *curNode  = effector->getParent();
    while (curNode != NULL) {
      LinkState* linkState = prevNode->getLinkState();
      vector<size_t> jointIds = linkState->getJointIdentifiers();
      assert(jointIds.size() == linkState->dof());
      ctx.popContext();
      vec3 direction = effectorPos - ctx.getCurrentOrigin();
      //cout << "dof: " << linkState->dof() << endl;
      //cout << "direction: " << direction << endl;
      for (size_t jointIdx = 0; jointIdx < linkState->dof(); jointIdx++) {
        size_t jointId = jointIds[jointIdx];
        //cout << "jointId: " << jointId << endl;
        vec3 rotAxis = linkState->getRotationAxis(jointIdx, ctx, desired.rows(3 * effectorId, 3 * effectorId + 2), effectorPos);
        //cout << "rotAxis: " << endl << rotAxis << endl;
        axes.rows(3 * jointId, 3 * jointId + 2) = rotAxis;
        vec3 jacobianEntry = cross(rotAxis, direction); 
        m(3 * effectorId,     jointId) = jacobianEntry[0];
        m(3 * effectorId + 1, jointId) = jacobianEntry[1];
        m(3 * effectorId + 2, jointId) = jacobianEntry[2];
      }
      prevNode = curNode;
      curNode  = curNode->getParent();
    }
  }

  //cout << "---" << endl;

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
  buf.resize(_effectors.size());
  for (vector<TreeNode*>::const_iterator it = _effectors.begin(); 
      it != _effectors.end();
      ++it) {
    size_t idx = (*it)->getIdentifier();
    buf[idx] = (*it)->getGlobalPosition(_rootPosition);
  }
  return buf;
}

vector<vec3>& LinkedTreeRobot::getInnerNodePositions(vector<vec3>& buf) const {
  buf.resize(_innerNodes.size());
  for (vector<TreeNode*>::const_iterator it = _innerNodes.begin(); 
      it != _innerNodes.end();
      ++it) {
    size_t idx = (*it)->getIdentifier();
    buf[idx] = (*it)->getGlobalPosition(_rootPosition);
  }
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

vec LinkedTreeRobot::computeDeltaThetas(const vec& desiredPositions, vec& axes) const {

  assert(_numEffectors * 3 == desiredPositions.n_elem);

  vec s(3 * _numEffectors);
  getEffectorPositions(s);

  vec e = clamp(desiredPositions - s, 0.5);

  mat J;
  computeJacobian(desiredPositions, J, axes);

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

void LinkedTreeRobot::updateThetas(const vec& deltas, const vec& axes) {
  assert(deltas.n_elem == _numJoints);
  assert(deltas.n_elem * 3 == axes.n_elem);
  _root->updateThetas(deltas, axes);
}

/** render links as lines for now */
void LinkedTreeRobot::renderRobot() const {
  Context ctx(_rootPosition);
  _root->renderTree(ctx);
}

void LinkedTreeRobot::solveIK(const vec& desired) {
  vec axes;
  vec deltaThetas = computeDeltaThetas(desired, axes);
  updateThetas(deltaThetas, axes);
}

}
}
