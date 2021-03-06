#ifndef _LINK_H_
#define _LINK_H_

#include <vector>
#include <armadillo>
#include <cassert>
#include <iostream>

#include "context.h"
#include "../util/buffer.h"

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

  /** set length */
  inline void setLength(const double);

  /** Given a context, compute the GLOBAL endpoint of this joint */
  inline arma::vec3 getEndpoint(const Context&) const;

  /** Degrees of freedom */
  virtual size_t dof() const = 0;

  /** Compute all jacobian entries */
  virtual void computeJacobianEntries(arma::mat&, /* The jacobian */
                                      arma::vec&, /* The state matrix */ 
                                      const size_t, /* Effector ID */
                                      const Context&, /* ctx */
                                      const arma::vec3&, /* desired effector pos */
                                      const arma::vec3&, /* current effector pos */
                                      const arma::vec3& /* current vector from origin TO effector */) const = 0;

  ///** For the i-th DOF, return the rotation axis given a context, desired
  // * position, and current effector position.
  // * DOF must be a valid DOF (indexed from zero). */
  //virtual arma::vec3 getRotationAxis(const size_t, 
  //                                   const Context&, 
  //                                   const arma::vec3&,
  //                                   const arma::vec3&) const = 0;

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

  /** Update the angles based on the global update vector. the rotation 
   * vector used is given as a reference */
  virtual void updateThetas(const arma::vec&, const arma::vec&) = 0;

  /** Set the angles based on the global update vector. the rotation 
   * vector used is given as a reference */
  virtual void setThetas(const arma::vec&, const arma::vec&) = 0;

  /** Get an orthonormal basis (u, v, n) for this frame */
  virtual void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const = 0;

  /** Mark a set point for this link. subsequent calls to reset will move the
   * link back to this exact state. calling mark overrides the previous set
   * point. note that mark has NOTHING to do with the animation buffer */
  virtual void mark() = 0;

  /** reset this link to the previously MARKED position. has nothing to
   * do with the animation buffer */
  virtual void reset() = 0;

                          /* Animation */

  /** Saves the current config to the animation buffer. does NOT affect
   * mark */
  virtual void save() = 0;

  /** restores this link to the current state pointed to by the animation
   * buffer pointer */
  virtual void restore() = 0;

  /** resets the animation buffer pointer to the beginning */
	virtual void resetPointer() = 0;

  /** clears the animation buffer and resets the pointer */
	virtual void resetBuffer() = 0;

	/** advance the animation buffer pointer */
	virtual void advance() = 0;

private:
  double _length;
};

inline double LinkState::getLength() const { return _length; }

inline void LinkState::setLength(const double _l) { _length = _l; }

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
  void computeJacobianEntries(arma::mat&,
                              arma::vec&,
                              const size_t,
                              const Context&,
                              const arma::vec3&,
                              const arma::vec3&,
                              const arma::vec3&) const;
  //arma::vec3 getRotationAxis(const size_t, 
  //                           const Context&, 
  //                           const arma::vec3&,
  //                           const arma::vec3&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void setThetas(const arma::vec&, const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

  void mark();
  void reset();

  void save();
  void restore();
	void resetPointer();
	void resetBuffer();
	void advance();

private:
  arma::vec3 axis; /* NORMALIZED axis of rotation */ 
  double angle; /* radians */
  arma::vec3 baselineDirection; /* When angle = 0, what NORMALIZED dir is this link pointing in? */

  arma::vec3 oldAxis;
  double oldAngle;
  arma::vec3 oldBaselineDirection;

  edu_berkeley_cs184::util::PlaybackBufferContainer3<arma::vec3, double, arma::vec3> _buf;

  size_t jointId;
};

class EulerBallAndSocketJoint : public LinkState {
public:
  EulerBallAndSocketJoint(const double, const arma::vec3&);
                  
  size_t dof() const;
  void computeJacobianEntries(arma::mat&,
                              arma::vec&,
                              const size_t,
                              const Context&,
                              const arma::vec3&,
                              const arma::vec3&,
                              const arma::vec3&) const;
  //arma::vec3 getRotationAxis(const size_t, 
  //                           const Context&, 
  //                           const arma::vec3&,
  //                           const arma::vec3&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void setThetas(const arma::vec&, const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

  void mark();
  void reset();

  void save();
  void restore();
	void resetPointer();
	void resetBuffer();
	void advance();

private:
  arma::vec3 baselineDirection;

  // Euler angles
  double thetax;
  double thetay;
  double thetaz;

  double oldThetax;
  double oldThetay;
  double oldThetaz;

  std::vector<size_t> jointIds;

};

class AxisBallAndSocketJoint : public LinkState {
public:
  AxisBallAndSocketJoint(const double, const arma::vec3&);

  size_t dof() const;
  void computeJacobianEntries(arma::mat&,
                              arma::vec&,
                              const size_t,
                              const Context&,
                              const arma::vec3&,
                              const arma::vec3&,
                              const arma::vec3&) const;
  //arma::vec3 getRotationAxis(const size_t, 
  //                           const Context&, 
  //                           const arma::vec3&,
  //                           const arma::vec3&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void setThetas(const arma::vec&, const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

  void mark();
  void reset();

  void save();
  void restore();
	void resetPointer();
	void resetBuffer();
	void advance();

private:
  arma::vec3 currentDirection; /** NORMALIZED */
  arma::mat44 currentTransform; /** local frame transform matrix */

  arma::vec3 oldDirection;
  arma::mat44 oldTransform;

  edu_berkeley_cs184::util::PlaybackBufferContainer2<arma::vec3, arma::mat44> _buf;

  size_t jointId;

};

class TranslationJoint : public LinkState {
public:
  TranslationJoint(const double, const arma::vec3&);

  size_t dof() const;
  void computeJacobianEntries(arma::mat&,
                              arma::vec&,
                              const size_t,
                              const Context&,
                              const arma::vec3&,
                              const arma::vec3&,
                              const arma::vec3&) const;
  //arma::vec3 getRotationAxis(const size_t, 
  //                           const Context&, 
  //                           const arma::vec3&,
  //                           const arma::vec3&) const; 
  arma::vec3 getRotatedDirection(const Context&) const; 
  std::vector<size_t> getJointIdentifiers() const;
  size_t assignJointIndicies(const size_t);
  Context& pushContext(Context&) const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void setThetas(const arma::vec&, const arma::vec&);
  void getBasis(const Context&, arma::vec3&, arma::vec3&, arma::vec3&) const;

  void mark();
  void reset();

  void save();
  void restore();
	void resetPointer();
	void resetBuffer();
	void advance();

private:
  arma::vec3 _fixedAxis; /** NORMALIZED */

  edu_berkeley_cs184::util::PlaybackBuffer<double> _buf;

  double oldLength;

  size_t jointId;
};


}
}

#endif /* _LINK_H_ */
