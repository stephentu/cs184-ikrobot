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

  _oldRootPosition = _rootPosition;
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
      //for (size_t jointIdx = 0; jointIdx < linkState->dof(); jointIdx++) {
      //  size_t jointId = jointIds[jointIdx];
      //  //cout << "jointId: " << jointId << endl;
      //  vec3 rotAxis = linkState->getRotationAxis(jointIdx, ctx, desired.rows(3 * effectorId, 3 * effectorId + 2), effectorPos);
      //  //cout << "rotAxis: " << endl << rotAxis << endl;
      //  axes.rows(3 * jointId, 3 * jointId + 2) = rotAxis;
      //  vec3 jacobianEntry = cross(rotAxis, direction); 
      //  m(3 * effectorId,     jointId) = jacobianEntry[0];
      //  m(3 * effectorId + 1, jointId) = jacobianEntry[1];
      //  m(3 * effectorId + 2, jointId) = jacobianEntry[2];
      //}
      
      linkState->computeJacobianEntries(m,
                                        axes,
                                        effectorId,
                                        ctx,
                                        desired.rows(3 * effectorId, 3 * effectorId + 2),
                                        effectorPos,
                                        direction);
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
      vec3 normDir = direction;
      normalize_vec3(normDir);
      for (size_t jointIdx = 0; jointIdx < linkState->dof(); jointIdx++) {
        size_t jointId = jointIds[jointIdx];
        vec3 rotAxis = axes.rows(3 * jointId, 3 * jointId + 2);
        vec3 dPdq; 
        if (vec3_equals(rotAxis, normDir)) // HACK: this is how we do constraints for translation joints
          dPdq = rotAxis;
        else
          dPdq = cross(rotAxis, direction);
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

static const double lambda = 1.6;
static const double lambda_squared = lambda * lambda;

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
  const double MAX_ITER = 100;
  size_t itersFinished = 0;
  const double TOL = 1e-3; 
  vec deltaThetas;
  do {
    vec axes;
    deltaThetas = computeDeltaThetas(desired, axes);
    updateThetas(deltaThetas, axes);
  } while (++itersFinished < MAX_ITER && norm(deltaThetas, 2) >= TOL);
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

vec& LinkedTreeRobot::getQDot(const vec& desired, vec& qdot, vec& axes) {
  vec P(3 * _numEffectors);
  getEffectorPositions(P); // position of the effectors currently

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
    
    mat JcT = trans(Jc);

    mat A = Jc * JcT;
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

    vec dinv_ut_b(d.n_elem);
    dinv_ut_b.zeros();

    for (size_t d_idx = 0; d_idx < d.n_elem; d_idx++) {
      double elem = d[d_idx];
      if (!double_equals(elem, 0.0))
        dinv_ut_b[d_idx] = 1.0 / elem * utb[d_idx];
    }

    gc = JcT * (V * dinv_ut_b); // Jc^T * lambda
  }

  qdot = ga + gc; // qdot
  return qdot;
}

#define SIMPLE_EULER 1

void LinkedTreeRobot::solveIKWithConstraints(const vec& desired) {

  assert(_numEffectors * 3 == desired.n_elem);

  vec thisUpdate;
  thisUpdate.zeros(_numJoints);

  const double MAX_ITER = 1000;
  size_t itersFinished = 0;
  const double TOL = 5e-4; // stopping condition for ||thisUpdate||

#if SIMPLE_EULER
  double h = 0.01; // initial h
#else
  double h = 0.2; // initial h
  const double EPS = 0.02; // upper bound for error produced per unit of t 
#endif

  bool redo = false;
  //cout << "----" << endl;
  do {

#if SIMPLE_EULER
    vec curAxes, f_tn_yn;
    getQDot(desired, f_tn_yn, curAxes); // f(t_n, y_n), fills curAxes
    thisUpdate = h * f_tn_yn;
    updateThetas(thisUpdate, curAxes); // push robot ahead
#else
    vec curAxes, f_tn_yn;
    getQDot(desired, f_tn_yn, curAxes); // f(t_n, y_n), fills curAxes

    vec alpha1 = h * f_tn_yn; // alpha1 = h f(t_n, y_n)

    // alpha2 = h/2 f(t_n, y_n) + h/2 f(t_n + h/2, y_n + h/2 f(t_n, y_n))

    vec temp_update = (h / 2.0) * f_tn_yn;
    updateThetas(temp_update, curAxes); // shift robot
      vec tempAxes, f_tn_mid;
      getQDot(desired, f_tn_mid, tempAxes); // f(t_n + h/2, y_n + h/2 f(t_n, y_n))
    updateThetas(-temp_update, curAxes); // unshift robot

    vec alpha2 = temp_update + (h / 2.0) * f_tn_mid;

    double r = norm(alpha1 - alpha2, 2);

    if (r < EPS) { // accept A2 with y_{n+1} = 2 A2 - A1
      thisUpdate = 2.0 * alpha2 - alpha1;
      updateThetas(thisUpdate, curAxes); // push robot ahead
      redo = false;
    } else 
      redo = true;

    h = 0.9 * EPS / r * h; // update step size
#endif

  } while (++itersFinished < MAX_ITER && (redo || norm(thisUpdate, 2) >= TOL));
  //cout << "Solution took " << itersFinished << " iterations" << endl
  //      << "h = " << h << endl;
  //cout << "----" << endl;
}

void LinkedTreeRobot::toggleConstraint(const size_t idx) {
  assert(0 <= idx && idx < _allNodes.size());
  //if (_allNodes[idx]->isLeafNode()) return;
  if (_allNodes[idx]->isFixed()) 
    _allNodes[idx]->setFixed(false);
  else {
    _allNodes[idx]->setFixed(true);
    _allNodes[idx]->setFixedPosition(_allNodes[idx]->getGlobalPosition(_rootPosition));
  }
}



}
}
