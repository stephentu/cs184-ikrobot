#include "link.h"
#include "../util/util.h"

#include <stdexcept>

#define TODO (throw runtime_error("TODO"))

using namespace std;
using namespace arma;
using namespace edu_berkeley_cs184::util;

namespace edu_berkeley_cs184 {
namespace robot {

LinkState::~LinkState() {}

RotationJoint::RotationJoint(const double _l,
                             const vec3& _axis,
                             const vec3& _baseline) :
  LinkState(_l), axis(_axis), angle(0.0), baselineDirection(_baseline) {

  // normalize axis and baseline just to be sure
  normalize_vec3(axis); 
  normalize_vec3(baselineDirection);

  // now assert that the dot product of axis and baselineDirection is zero
  assert( double_equals(dot(axis, baselineDirection), 0) );

  mark();
}

size_t RotationJoint::dof() const { return 1; }

void RotationJoint::computeJacobianEntries(mat& jacobian,
                                           vec& state,
                                           const size_t effectorId,
                                           const Context& ctx,
                                           const vec3& desired,
                                           const vec3& effector,
                                           const vec3& direction) const {
  vec3 rotAxis = ctx.getVectorInContext(axis);
  state.rows(3 * jointId, 3 * jointId + 2) = rotAxis;
  vec3 jacobianEntry = cross(rotAxis, direction); 
  jacobian(3 * effectorId,     jointId) = jacobianEntry[0];
  jacobian(3 * effectorId + 1, jointId) = jacobianEntry[1];
  jacobian(3 * effectorId + 2, jointId) = jacobianEntry[2];
}

vec3 RotationJoint::getRotatedDirection(const Context& ctx) const {
  //cout << "axis: " << ctx.getVectorInContext(axis) << endl;
  return rotate_expmap(
    ctx.getVectorInContext(baselineDirection), 
    ctx.getVectorInContext(axis), 
    angle);
}

vector<size_t> RotationJoint::getJointIdentifiers() const {
  vector<size_t> v;
  v.push_back(jointId);
  return v;
}

size_t RotationJoint::assignJointIndicies(const size_t startIdx) {
  jointId = startIdx;
  return 1;
}

Context& RotationJoint::pushContext(Context& ctx) const {
  ctx.pushContext(getEndpoint(ctx), axis, angle);
  return ctx;
}

void RotationJoint::updateThetas(const vec& deltas, const vec& axes) {
  angle += deltas(jointId);
}

void RotationJoint::setThetas(const vec& deltas, const vec& axes) {
  angle = deltas(jointId);
}

void RotationJoint::getBasis(const Context& ctx, vec3& u, vec3& v, vec3& n) const {
  n = getRotatedDirection(ctx);
  u = ctx.getVectorInContext(axis); 
  v = cross(n, u);
  normalize_vec3(v);
}

void RotationJoint::toggleConstraint() { TODO; }

void RotationJoint::getConstraintInfo(
      arma::vec& constraints,
      arma::mat& constraintGradient,
      const Context& ctx,
      const arma::vec& jointIds,
      const arma::vec& rotationAxes
    ) {
  TODO;
}

size_t RotationJoint::numConstraints() const {
  return 0;
}

void RotationJoint::mark() {
  oldAxis              = axis;
  oldAngle             = angle;
  oldBaselineDirection = baselineDirection;
}

void RotationJoint::reset() {
  axis              = oldAxis;
  angle             = oldAngle;
  baselineDirection = oldBaselineDirection;
}

void RotationJoint::save() { 
  _buf.append(axis, angle, baselineDirection);
}

void RotationJoint::restore() {
  Tuple3<vec3, double, vec3> values = _buf.getCurrent();
  axis              = values._1;
  angle             = values._2;
  baselineDirection = values._3;
}

void RotationJoint::resetPointer() { _buf.resetPointer(); }
void RotationJoint::resetBuffer() { _buf.resetBuffer(); }
void RotationJoint::advance() { _buf.advance(); }

AxisBallAndSocketJoint::AxisBallAndSocketJoint(const double _l, const vec3& startingDir) 
  : LinkState(_l), currentDirection(startingDir) {
  normalize_vec3(currentDirection);
  currentTransform.eye();
  mark();
}

size_t AxisBallAndSocketJoint::dof() const { return 1; }

void AxisBallAndSocketJoint::computeJacobianEntries(mat& jacobian,
                                                    vec& state,
                                                    const size_t effectorId,
                                                    const Context& ctx,
                                                    const vec3& desired,
                                                    const vec3& effector,
                                                    const vec3& direction) const {
  vec3 rotAxis;
  if (vec3_equals(desired, effector))
    rotAxis = makeVec3(1, 0, 0);
  else {
    rotAxis = cross(effector - ctx.getCurrentOrigin(),
                    desired - ctx.getCurrentOrigin());
    rotAxis = normalize_vec3(rotAxis);
  }

  state.rows(3 * jointId, 3 * jointId + 2) = rotAxis;
  vec3 jacobianEntry = cross(rotAxis, direction); 
  jacobian(3 * effectorId,     jointId) = jacobianEntry[0];
  jacobian(3 * effectorId + 1, jointId) = jacobianEntry[1];
  jacobian(3 * effectorId + 2, jointId) = jacobianEntry[2];
}

vec3 AxisBallAndSocketJoint::getRotatedDirection(const Context& ctx) const {
  return ctx.getVectorInContext(currentDirection);
}

vector<size_t> AxisBallAndSocketJoint::getJointIdentifiers() const {
  vector<size_t> v;
  v.push_back(jointId);
  return v;
}

size_t AxisBallAndSocketJoint::assignJointIndicies(const size_t startIdx) {
  jointId = startIdx;
  return 1;
}

Context& AxisBallAndSocketJoint::pushContext(Context& ctx) const {
  ctx.pushContext(getEndpoint(ctx), currentTransform);
  return ctx;
}

void AxisBallAndSocketJoint::updateThetas(const vec& deltas, const vec& axes) {
  const double dt = deltas(jointId);
  const vec3 axis = axes.rows(3 * jointId, 3 * jointId + 2);

  currentDirection = rotate_expmap(currentDirection, axis, dt); 
  currentTransform = rotation_matrix(axis, dt) * currentTransform; 
}

void AxisBallAndSocketJoint::setThetas(const vec& deltas, const vec& axes) {
  throw runtime_error("setThetas unsupported for AxisBallAndSocketJoint");
}

void AxisBallAndSocketJoint::getBasis(const Context& ctx, vec3& u, vec3& v, vec3& n) const {
  n = ctx.getVectorInContext(currentDirection);
  if (vec3_equals(makeVec3(0, 1, 0), n) || vec3_equals(makeVec3(0, -1, 0), n))
    u = cross(makeVec3(1, 0, 0), n);
  else
    u = cross(makeVec3(0, 1, 0), n);
  v = cross(n, u);
  normalize_vec3(u);
  normalize_vec3(v);
  normalize_vec3(n);
}

void AxisBallAndSocketJoint::toggleConstraint() { TODO; }

void AxisBallAndSocketJoint::getConstraintInfo(
      arma::vec& constraints,
      arma::mat& constraintGradient,
      const Context& ctx,
      const arma::vec& jointIds,
      const arma::vec& rotationAxes
    ) {
  TODO;
}

size_t AxisBallAndSocketJoint::numConstraints() const {
  return 0;
}

void AxisBallAndSocketJoint::mark() {
  oldDirection = currentDirection;
  oldTransform = currentTransform;
}

void AxisBallAndSocketJoint::reset() {
  currentDirection = oldDirection;
  currentTransform = oldTransform;
}

void AxisBallAndSocketJoint::save() { 
  _buf.append(currentDirection, currentTransform);
}

void AxisBallAndSocketJoint::restore() { 
  pair<vec3, mat44> values = _buf.getCurrent();
  currentDirection = values.first;
  currentTransform = values.second;
}

void AxisBallAndSocketJoint::resetPointer() { _buf.resetPointer(); }
void AxisBallAndSocketJoint::resetBuffer() { _buf.resetBuffer(); }
void AxisBallAndSocketJoint::advance() { _buf.advance(); }

TranslationJoint::TranslationJoint(const double initLen, const vec3& dir) 
  : LinkState(initLen), _fixedAxis(dir) {
  normalize_vec3(_fixedAxis);
  mark();
}

size_t TranslationJoint::dof() const { return 1; }

void TranslationJoint::computeJacobianEntries(mat& jacobian,
                                              vec& state,
                                              const size_t effectorId,
                                              const Context& ctx,
                                              const vec3& desired,
                                              const vec3& effector,
                                              const vec3& direction) const {
  jacobian(3 * effectorId,     jointId) = _fixedAxis[0];
  jacobian(3 * effectorId + 1, jointId) = _fixedAxis[1];
  jacobian(3 * effectorId + 2, jointId) = _fixedAxis[2];
}

vec3 TranslationJoint::getRotatedDirection(const Context& ctx) const {
  return ctx.getVectorInContext(_fixedAxis);
}

vector<size_t> TranslationJoint::getJointIdentifiers() const {
  vector<size_t> v;
  v.push_back(jointId);
  return v;
}

size_t TranslationJoint::assignJointIndicies(const size_t startIdx) {
  jointId = startIdx;
  return 1;
}

Context& TranslationJoint::pushContext(Context& ctx) const {
  ctx.pushContext(getEndpoint(ctx));
  return ctx;
}

void TranslationJoint::updateThetas(const vec& deltas, const vec& axes) {
  const double dt = deltas(jointId);
  setLength(getLength() + dt);
}

void TranslationJoint::setThetas(const vec& deltas, const vec& axes) {
  const double dt = deltas(jointId);
  setLength(dt);
}

void TranslationJoint::getBasis(const Context& ctx, vec3& u, vec3& v, vec3& n) const {
  n = ctx.getVectorInContext(_fixedAxis);
  if (vec3_equals(makeVec3(0, 1, 0), n) || vec3_equals(makeVec3(0, -1, 0), n))
    u = cross(makeVec3(1, 0, 0), n);
  else
    u = cross(makeVec3(0, 1, 0), n);
  v = cross(n, u);
  normalize_vec3(u);
  normalize_vec3(v);
  normalize_vec3(n);
}

void TranslationJoint::toggleConstraint() {

}

void TranslationJoint::getConstraintInfo(
      arma::vec& constraints,
      arma::mat& constraintGradient,
      const Context& ctx,
      const arma::vec& jointIds,
      const arma::vec& rotationAxes
    ) {
  constraints.set_size(1);
  double x = getLength() - 1.0;
  constraints[0] = x * x;
  constraintGradient.zeros(1, jointIds.n_elem);
  constraintGradient[jointId] = 2.0 * x;
}

size_t TranslationJoint::numConstraints() const {
  return 1;
}

void TranslationJoint::mark() {
  oldLength = getLength();
}

void TranslationJoint::reset() {
  setLength(oldLength);
}

void TranslationJoint::save() { 
  _buf.append(getLength());
}

void TranslationJoint::restore() { 
  setLength(_buf.getCurrent());
}

void TranslationJoint::resetPointer() { _buf.resetPointer(); }
void TranslationJoint::resetBuffer() { _buf.resetBuffer(); }
void TranslationJoint::advance() { _buf.advance(); }

}
}
