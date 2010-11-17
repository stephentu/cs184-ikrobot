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

  _root->gatherNodes(_allNodes);

  assert(nJoints == _numJoints);
  assert(nLeaves == _numEffectors);
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
    size_t effectorId = effector->getLeafIdentifier();

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

mat& LinkedTreeRobot::computeConstraintJacobian(const vec& desired,
                                                const vec& axes,
                                                mat& m) const {

  vector<TreeNode*> constraintNodes;
  for (vector<TreeNode*>::const_iterator it = _allNodes.begin();
      it != _allNodes.end();
      ++it)
    if ((*it)->isFixed())
      constraintNodes.push_back(*it);

  m.zeros(3 * constraintNodes.size(), _numJoints);

  for (size_t idx = 0; idx < constraintNodes.size(); idx++) {
    TreeNode *node = constraintNodes[idx];

    Context ctx(_rootPosition);
    node->getContextForNode(ctx);

    vec3 currentPos = ctx.getCurrentOrigin();
    vec3 constraintError = node->getFixedPosition() - currentPos;
    vec3 dCdP = -2.0 * constraintError;

    TreeNode *curNode = node;
    while (!curNode->isRootNode()) {
      //cout << "currentNode: " << curNode->getIdentifier() << endl;
      LinkState* linkState = curNode->getLinkState();
      vector<size_t> jointIds = linkState->getJointIdentifiers();
      assert(jointIds.size() == linkState->dof());
      ctx.popContext();
      vec3 direction = currentPos - ctx.getCurrentOrigin();
      for (size_t jointIdx = 0; jointIdx < linkState->dof(); jointIdx++) {
        size_t jointId = jointIds[jointIdx];
        vec3 rotAxis = axes.rows(3 * jointId, 3 * jointId + 2);
        vec3 dPdq = cross(rotAxis, direction);
        vec3 jacEntry = dCdP % dPdq;
        m.submat(3 * idx, jointId, 3 * idx + 2, jointId) = jacEntry; 
      }
      curNode = curNode->getParent();
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
  buf.resize(_effectors.size());
  for (vector<TreeNode*>::const_iterator it = _effectors.begin(); 
      it != _effectors.end();
      ++it) {
    size_t idx = (*it)->getLeafIdentifier();
    buf[idx] = (*it)->getGlobalPosition(_rootPosition);
  }
  return buf;
}

vector<vec3>& LinkedTreeRobot::getNodePositions(vector<vec3>& buf) const {
  buf.resize(_allNodes.size());
  for (vector<TreeNode*>::const_iterator it = _allNodes.begin(); 
      it != _allNodes.end();
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

static inline bool is_all_zeros(const mat& m) {
  for (mat::const_iterator it = m.begin();
      it != m.end(); ++it)
    if (!double_equals(*it, 0.0))
      return false;
  return true;
}

static inline mat truncate(const mat& m, const double threshold) {
  mat mtrunc;
  mtrunc.zeros(m.n_rows, m.n_cols);
  for (size_t i = 0; i < m.n_rows; i++)
    for (size_t j = 0; j < m.n_cols; j++)
      if (m(i, j) >= threshold)
        mtrunc(i, j) = m(i, j);
  return mtrunc;
}

void LinkedTreeRobot::solveIKWithConstraints(const vec& desired) {
  vec axes;
  assert(_numEffectors * 3 == desired.n_elem);

  vec lastUpdate, thisUpdate;
  lastUpdate.zeros(_numJoints);
  thisUpdate.zeros(_numJoints);

  const double TOL = 1e-6;
  size_t itersFinished = 0;

  do {
    lastUpdate = thisUpdate;

    vec P(3 * _numEffectors);
    getEffectorPositions(P);

    vec F = clamp(desired - P, 0.5);

    mat J;
    computeJacobian(desired, J, axes);

    mat Jc;
    //cout << "++ Going to compute constraint jacobian ++" << endl;
    computeConstraintJacobian(desired, axes, Jc);
    //cout << "++ Finished computed constraint jacobian ++" << endl;
    //cout << "Jc:" << endl << Jc << endl;

    vec ga = trans(J) * F;
    //cout << "ga:" << endl << ga << endl;

    vec gc;
    if (Jc.n_elem == 0 || is_all_zeros(Jc)) // no constraints, or all constraints OK
      gc.zeros(ga.n_elem);
    else {
      vec b = -Jc * ga;
      //cout << "b:" << endl << b << endl;

      mat A = Jc * trans(Jc);
      //cout << "A:" << endl << A << endl;

      vec d;
      mat U, V;
      if (!svd(U, d, V, A))
        cerr << "could not compute SVD" << endl;

      //cout << "U:" << endl << U << endl;
      //cout << "V:" << endl << V << endl;
      //cout << "d:" << endl << d << endl;

      //vec lambda = V * inv(D) * trans(U) * b; // too naive

      vec utb = trans(U) * b;
      vec dinvcols;
      dinvcols.zeros(d.n_elem);
      for (size_t d_idx = 0; d_idx < d.n_elem; d_idx++)
        if (!double_equals(d(d_idx), 0.0))
          dinvcols(d_idx) = 1.0 / d(d_idx);

      vec dinv_ut_b = diagmat(dinvcols) * utb; 

      gc = trans(Jc) * (V * dinv_ut_b); // Jc^T * lambda
    }

    vec qdot = ga + gc;

    thisUpdate = 0.01 * qdot; // forward euler

    updateThetas(thisUpdate, axes);

    //if ((++itersFinished % 100) == 0) 
    //  cout << "Just finished " << itersFinished << " iterations" << endl;
  } while (norm(thisUpdate - lastUpdate, 2) >= TOL);
}

void LinkedTreeRobot::toggleConstraint(const size_t idx) {
  assert(0 <= idx && idx < _allNodes.size());
  if (_allNodes[idx]->isFixed()) 
    _allNodes[idx]->setFixed(false);
  else {
    _allNodes[idx]->setFixed(true);
    _allNodes[idx]->setFixedPosition(_allNodes[idx]->getGlobalPosition(_rootPosition));
  }
}



}
}
