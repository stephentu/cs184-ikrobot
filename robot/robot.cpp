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

static const double lambda = 1.5;
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

vec& LinkedTreeRobot::computeGradL(
    vec& gradL, /** gradL */
    vec& rotAxes, /** state vector */
    const vec& desired, /* desired positions */
    const vec& S, /* S */
    const vec& lambda, /* lagrange multipliers */
    const vector<LinkState*>& linkStates, /* all link states */
    const size_t numConstraints
  ) {

  assert(desired.n_elem == 3 * _numEffectors);
  assert(S.n_elem == numConstraints);
  assert(lambda.n_elem == numConstraints);

  vec q(3 * _numEffectors);
  getEffectorPositions(q);

  vec e = clamp(desired - q, 0.5);

  //cout << "desired: " << endl << desired << endl;
  //cout << "q: " << endl << q << endl;
  //cout << "S: " << endl << S << endl;
  //cout << "lambda: " << endl << lambda << endl;

  mat J; // d(effectors)/dq
  computeJacobian(desired, J, rotAxes);
  //cout << "J: " << endl << J << endl;

  vec dLdq, dLds, dLdlambda;
  //dLds = -2.0 * (lambda % S);
  dLds.zeros(numConstraints);
  dLdlambda.zeros(numConstraints);

  vec gradf;
  gradf.zeros(_numJoints);
  for (size_t i = 0; i < _numJoints; i++) {
    //cout << "e % J.col(i)" << endl << (e % J.col(i)) << endl;
    double x = -as_scalar(sum(e % J.col(i)));
    //cout << "x = " << x << endl;
    gradf(i) = x;
  }

  //cout << "grad F = " << endl << gradf << endl;

  dLdq = gradf; // assumes grad(C(q)) == 0
  size_t constraintIdx = 0;
  for (vector<LinkState*>::const_iterator it = linkStates.begin();
      it != linkStates.end();
      ++it) {
    LinkState *ls = *it;

    Context ctx(_rootPosition);
    vec constraints;
    mat constraintsGrad;
    vec jointIds, rotationAxes; // TODO: fill in 

    jointIds.set_size(1);
    assert(ls->dof() == 1);
    jointIds[0] = ls->getJointIdentifiers()[0];

    size_t thisNumConstraints = ls->numConstraints();
    ls->getConstraintInfo(constraints, constraintsGrad, ctx, jointIds, rotationAxes);

    //dLdlambda.rows(constraintIdx, constraintIdx + thisNumConstraints - 1) =
    //  - 1.0 * ( (S.rows(constraintIdx, constraintIdx + thisNumConstraints - 1) %
    //             S.rows(constraintIdx, constraintIdx + thisNumConstraints - 1)) +
    //            constraints   
    //          );

    dLdlambda.rows(constraintIdx, constraintIdx + thisNumConstraints - 1) =
      - constraints;

    //cout << "constraints: " << endl << constraints << endl;
    //cout << "constraintsGrad: " << endl << constraintsGrad << endl;

    for (size_t constraintId = 0; constraintId < ls->numConstraints(); constraintId++) {
      vec gradC = trans(constraintsGrad.row(constraintId));
      double lambdai = dLdlambda(constraintIdx + constraintId);
      //cout << "lambda_i: " << lambdai << endl;
      //cout << "gradC: " << endl << gradC << endl;
      //cout << "prod: " << (lambdai * gradC) << endl;
      dLdq = dLdq - lambdai * gradC;
    }

    constraintIdx += thisNumConstraints;
  }

  gradL = join_cols(join_cols(dLdq, dLds), dLdlambda);
  assert(gradL.n_elem ==  _numJoints + 2 * numConstraints);

  //cout << "grad L = " << endl << gradL << endl;

  return gradL;
}

void LinkedTreeRobot::solveIKGradientDescent(const vec& desired) {
  // use lagrange function:
  // L(q, s, lambda) = f(q) - lambda * (c(q) + s^2)
  // find min L(q, s, lambda)
  // f(q) = || desired - currentPositions(q) ||^2

  vector<LinkState*> linkStates;
  _root->gatherLinkStates(linkStates);

  size_t numConstraints = 0;
  for (vector<LinkState*>::iterator it = linkStates.begin();
      it != linkStates.end();
      ++it)
    numConstraints += (*it)->numConstraints();

  vec S, lambda, deltaS, deltaLambda, deltaThetas;
  S.ones(numConstraints);
  lambda.zeros(numConstraints);
  deltaThetas.zeros(_numJoints);

  const double MAX_ITER = 100000000; 
  size_t itersFinished = 0;
  const double TOL = 1e-3;
  const double STEP_SIZE = 0.0001;

  vec axes, gradL;
  do {
    computeGradL(gradL, axes, desired, S, lambda, linkStates, numConstraints);
    deltaThetas = - STEP_SIZE * gradL.rows(0, _numJoints - 1);
    deltaS = - STEP_SIZE * gradL.rows(_numJoints, _numJoints + numConstraints - 1);
    deltaLambda = - STEP_SIZE * gradL.rows(_numJoints + numConstraints, gradL.n_elem - 1);
    updateThetas(deltaThetas, axes);
    S += deltaS;
    lambda += deltaLambda;
  } while (++itersFinished < MAX_ITER && norm(gradL, 2) >= TOL);
  cout << "Finished in " << itersFinished << " iterations" << endl;
}

void LinkedTreeRobot::toggleConstraint(const size_t idx) {
  assert(0 <= idx && idx < _allNodes.size());

}



}
}
