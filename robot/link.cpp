#include "link.h"

#include "../util/util.h"

#include <stdexcept>

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
}

size_t RotationJoint::dof() const { return 1; }

vec3 RotationJoint::getRotationAxis(const size_t dof, const Context& ctx) const {
  assert(dof == 0);
  return ctx.getVectorInContext(axis);
}

vec3 RotationJoint::getRotatedDirection(const Context& ctx) const {
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

void RotationJoint::updateThetas(const vec& deltas) {
  angle += deltas(jointId);
}

void RotationJoint::getBasis(const Context& ctx, vec3& u, vec3& v, vec3& n) const {
  n = getRotatedDirection(ctx);
  u = getRotationAxis(0, ctx);
  v = cross(n, u);
}

BallAndSocketJoint::BallAndSocketJoint(const double _l, const vec3& _baseline) :
  LinkState(_l), baselineDirection(_baseline), thetax(0.0), thetay(0.0), thetaz(0.0) {

  normalize_vec3(baselineDirection);

}

size_t BallAndSocketJoint::dof() const { return 3; }

vec3 BallAndSocketJoint::getRotationAxis(const size_t dof, const Context& ctx) const {
  assert(0 <= dof && dof < 3);
  switch (dof) {
  case 0:
    return makeVec3(1, 0, 0);
  case 1:
    return makeVec3(0, 1, 0);
  case 2:
    return makeVec3(0, 0, 1);
  default:
    throw runtime_error("illegal DOF");
  }
}

vec3 BallAndSocketJoint::getRotatedDirection(const Context& ctx) const {
  return rotate_euler(ctx.getVectorInContext(baselineDirection), thetax, thetay, thetaz);
}

vector<size_t> BallAndSocketJoint::getJointIdentifiers() const { return jointIds; }

size_t BallAndSocketJoint::assignJointIndicies(const size_t startIdx) {
  jointIds.push_back(startIdx);
  jointIds.push_back(startIdx + 1);
  jointIds.push_back(startIdx + 2);
  return 3;
}

Context& BallAndSocketJoint::pushContext(Context& ctx) const {
  ctx.pushContext(getEndpoint(ctx), thetax, thetay, thetaz);
  return ctx;
}

void BallAndSocketJoint::updateThetas(const vec& deltas) {
  thetax += deltas(jointIds[0]);
  thetay += deltas(jointIds[1]);
  thetaz += deltas(jointIds[2]);
}

void BallAndSocketJoint::getBasis(const Context& ctx, vec3& u, vec3& v, vec3& n) const {
  n = getRotatedDirection(ctx);
  u = cross(makeVec3(0, 1, 0), n); // TODO: this might be degenerate sometimes
  v = cross(n, u);
}

}
}
