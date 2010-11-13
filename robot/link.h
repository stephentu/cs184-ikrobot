#ifndef _LINK_H_
#define _LINK_H_

#include <vector>
#include <armadillo>
#include <cassert>
#include <iostream>

#include "context.h"

namespace edu_berkeley_cs184 {
namespace robot {

class LinkState {
public:
  LinkState(const double _l) : _length(_l) {
    assert(_length > 0.0);
  }
  virtual ~LinkState();

  /** length of joint */
  inline double getLength() const;

  /** Given a context, compute the GLOBAL endpoint of this joint */
  inline arma::vec3 getEndpoint(const Context&) const;

  /** Degrees of freedom */
  virtual size_t dof() const = 0;

  /** For the i-th DOF, return the rotation axis given a context
   * DOF must be a valid DOF (indexed from zero) */
  virtual arma::vec3 getRotationAxis(const size_t, const Context&) const = 0;

  /** what NORMALIZED direction is the link pointing in, given a context and applying
   * this link's rotations?  */
  virtual arma::vec3 getRotatedDirection(const Context&) const = 0;

  /** get a joint ID for each DOF */
  virtual std::vector<size_t> getJointIdentifiers() const = 0;

  /** assign joint indicies, starting with the input, returning the number of
   * joint indicies assigned ( = DOF ) */
  virtual size_t assignJointIndicies(const size_t) = 0;

  /** push this joint's configuration onto the current context */
  virtual Context& pushContext(Context&) const = 0;

  /** Update the angles based on the global update vector */
  virtual void updateThetas(const arma::vec&) = 0;

  /** Get an orthonormal basis (u, v, n) for this frame */
  virtual void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const = 0;

private:
  double _length;
};

inline double LinkState::getLength() const { return _length; }

inline arma::vec3 LinkState::getEndpoint(const Context& ctx) const {
  // step 1: rotate baselineDirection by angle degrees along axis vector (in
  // right hand rule sense)
  arma::vec3 rotated = getRotatedDirection(ctx);

  //std::cout << "LinkState::getEndpoint: rotated: " << std::endl << rotated << std::endl;

  // step 2: endpoint = jointLoc + length * rotated
  return ctx.getCurrentOrigin() + getLength() * rotated;
}

class RotationJoint: public LinkState {
public:
  /** new link state with length _l, GLOBAL rotation axis _axis (w.r.t theta =
   * 0 for all links), and GLOBAL theta = 0 direction (w.r.t all theta = 0) */
  RotationJoint(const double, const arma::vec3&, const arma::vec3&);
  
  size_t dof() const;
  arma::vec3 getRotationAxis(const size_t, const Context&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

private:
  arma::vec3 axis; /* NORMALIZED axis of rotation */ 
  double angle; /* radians */
  arma::vec3 baselineDirection; /* When angle = 0, what NORMALIZED dir is this link pointing in? */

  size_t jointId;
};

class BallAndSocketJoint : public LinkState {
public:
  BallAndSocketJoint(const double, const arma::vec3&);
                  
  size_t dof() const;
  arma::vec3 getRotationAxis(const size_t, const Context&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

private:
  arma::vec3 baselineDirection;

  // Euler angles
  double thetax;
  double thetay;
  double thetaz;

  std::vector<size_t> jointIds;

};


}
}

#endif /* _LINK_H_ */
