#include "context.h"

#include "../util/util.h"

using namespace std;
using namespace arma;
using namespace edu_berkeley_cs184::util;

namespace edu_berkeley_cs184 {
namespace robot {

Context::Context() {
  _positions.push(makeVec3(0, 0, 0));
  mat44 I4;
  I4.eye();
  _transforms.push(I4);
}

Context::Context(const vec3& pos) {
  _positions.push(pos);
  mat44 I4;
  I4.eye();
  _transforms.push(I4);
}

Context::Context(const vec3& pos,
                 const mat44& trfm) {
  _positions.push(pos);
  _transforms.push(trfm);
}

/** Push a new context of a new position, with a new rotation vector + the
 * amount of rotation (in radians) */
void Context::pushContext(const vec3& newPosition,
                          const vec3& newRotVector,
                          const double radians) {
  _positions.push(newPosition);
  _transforms.push(rotation_matrix(newRotVector, radians) * _transforms.top());
}

/** Push a new context of a new position + euler angles */
void Context::pushContext(const vec3& newPosition,
                          const double thx,
                          const double thy,
                          const double thz) {
  _positions.push(newPosition);
  _transforms.push(rotation_matrix(thx, thy, thz) * _transforms.top());
}

void Context::popContext() {
  _positions.pop();
  _transforms.pop();
}


}
}

