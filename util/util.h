#ifndef _UTIL_H_
#define _UTIL_H_

#include <armadillo>
#include <cmath>
#include <limits>
#include <cassert>

namespace edu_berkeley_cs184 {
namespace util {

inline arma::vec3& normalize_vec3(arma::vec3& v) {
  double x1 = v[0];
  double x2 = v[1];
  double x3 = v[2];
  double l = sqrt( x1 * x1 + x2 * x2 + x3 * x3 );
  v[0] = x1 / l;
  v[1] = x2 / l;
  v[2] = x3 / l;
  return v;
}

inline bool double_equals(const double x1, const double x2) {
  return fabs(x1 - x2) < 1000.0 * std::numeric_limits<double>::epsilon();
}

inline bool vec3_equals(const arma::vec3& a, const arma::vec3& b) {
  return double_equals(a[0], b[0]) &&
         double_equals(a[1], b[1]) &&
         double_equals(a[2], b[2]);
}


/**
 * axis must be normalized. theta must be in radians
 */
inline arma::vec3 rotate_expmap(const arma::vec3& v, 
                                const arma::vec3& axis,
                                const double theta) {

  // use the rodrigues formula to compute rotation. requires axis to be
  // normalized. 
  assert( double_equals(arma::norm(axis, 2), 1.0) );

  // v_rot = (v cos(theta)) + ((z x v) sin(theta)) + z(z dot v)(1 - cos(theta))
  // where z = axis

  const double cos_th = cos(theta);
  const double sin_th = sin(theta);

  return (v * cos_th) + ((cross(axis, v) * sin_th)) + axis * dot(axis, v) * (1.0 - cos_th);
}

inline arma::vec3 makeVec3(const double x1, 
                           const double x2, 
                           const double x3) {
  arma::vec3 v;
  v[0] = x1;
  v[1] = x2;
  v[2] = x3;
  return v;
}

inline arma::rowvec3 makeRowVec3(const double x1,
                                 const double x2,
                                 const double x3) {
  arma::rowvec3 v;
  v[0] = x1;
  v[1] = x2;
  v[2] = x3;
  return v;
}


inline arma::mat44 rotation_matrix(const arma::vec3& z, // axis to rotate about - assumes is already normal
                                   const double theta /* angle in radians */) {
  arma::mat33 zcross;
  zcross.row(0) = makeRowVec3(0, -z[2], z[1]);
  zcross.row(1) = makeRowVec3(z[2], 0, -z[0]);
  zcross.row(2) = makeRowVec3(-z[1], z[0], 0);

  arma::mat33 I3;
  I3.eye();

  const double cos_th = cos(theta);
  const double sin_th = sin(theta);

  arma::mat33 R = (I3 * cos_th + zcross * sin_th + (1.0 - cos_th) * z * trans(z)); 

  arma::mat44 R4;
  R4.eye();
  R4.submat(arma::span(0, 2), arma::span(0, 2)) = R;

  return R4;
}

/** Euler angles in radians. Treated as RxRyRz */
inline arma::mat44 rotation_matrix(const double thx,
                                   const double thy,
                                   const double thz) {

  const double cx = cos(thx);
  const double cy = cos(thy);
  const double cz = cos(thz);

  const double sx = sin(thx);
  const double sy = sin(thy);
  const double sz = sin(thz);

  arma::mat33 R;
  R.row(0) = makeRowVec3(cy * cz, -cy * sz, sy);
  R.row(1) = makeRowVec3(cz * sx * sy + cx * sz, cx * cz - sx * sy * sz, -cy * sx);
  R.row(2) = makeRowVec3(-cx * cz * sy + sx * sz, cz * sx + cx * sy * sz, cx * cy);

  arma::mat44 R4;
  R4.eye();
  R4.submat(arma::span(0, 2), arma::span(0, 2)) = R;

  return R4;
}

inline arma::vec3 rotate_euler(const arma::vec3& v,
                               const double thx,
                               const double thy,
                               const double thz) {
  arma::mat44 R = rotation_matrix(thx, thy, thz);
  return R.submat(arma::span(0, 2), arma::span(0, 2)) * v;
}


}
}

#endif /* _UTIL_H_ */
