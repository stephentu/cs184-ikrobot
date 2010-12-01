#ifndef _ROBOT_H_
#define _ROBOT_H_

#include <armadillo>
#include <iostream>

#include "tree.h"

namespace edu_berkeley_cs184 {
namespace robot {

enum SolnType {
  PINV,
  DLS,
  NUM_TYPES
};

class LinkedTreeRobot {
public:
  LinkedTreeRobot(const arma::vec3&, TreeNode*);
  virtual ~LinkedTreeRobot();

  inline size_t getNumJoints() const;
  inline size_t getNumEffectors() const;

  /** Get the effector positions as one large vector */
  arma::vec& getEffectorPositions(arma::vec&) const;

  /** Get effector positions as an stl vector of arma vec3s */
  std::vector<arma::vec3>& getEffectorPositions(std::vector<arma::vec3>&) const;

  /** Get each of the inner node positions as one large vector */
  std::vector<arma::vec3>& getNodePositions(std::vector<arma::vec3>&) const;

  arma::mat& computeJacobian(const arma::vec&, arma::mat&, arma::vec&) const;
  arma::mat& computeConstraintJacobian(const arma::vec&, const arma::vec&, arma::mat&) const;
  arma::vec computeDeltaThetas(const arma::vec&, arma::vec&) const;
  void updateThetas(const arma::vec&, const arma::vec&);

  /** Main IK method- given an input of desired positions, update the robot to
   * try to match this input */
  void solveIK(const arma::vec&);

  /** Main IK with constraints method. same input as solveIK. uses completely
   * different algorithm (iterative jacobian transpose with lagrange
   * multipliers). the method type is ignored when using this method */
  void solveIKWithConstraints(const arma::vec&);

  /** draw the robot in openGL */
  void renderRobot() const;

  inline SolnType getMethod() const;

  inline void setMethod(const SolnType);

  inline void setRootPosition(const arma::vec3&);

  void toggleConstraint(const size_t);

private:
  arma::vec& getQDot(const arma::vec&, arma::vec&, arma::vec&);

  arma::vec3 _rootPosition;
  TreeNode* _root;

  // computed from _root
  size_t _numJoints;
  size_t _numEffectors;

  // effectors, from traversing the root
  std::vector<TreeNode*> _effectors;

  // inner nodes, from traversing the root. this + effectors gives us the
  // entire tree from root
  std::vector<TreeNode*> _innerNodes;

  // all nodes from root
  std::vector<TreeNode*> _allNodes;

  // method to use when solving for delta thetas
  SolnType _method;
};

inline size_t LinkedTreeRobot::getNumJoints() const {
  return _numJoints;
}

inline size_t LinkedTreeRobot::getNumEffectors() const {
  return _numEffectors;
}

inline SolnType LinkedTreeRobot::getMethod() const { return _method; }
inline void LinkedTreeRobot::setMethod(const SolnType _t) { 
  std::cout << "Setting method to: " << _t << std::endl;
  _method = _t; 
}

inline void LinkedTreeRobot::setRootPosition(const arma::vec3& _p) { 
  //std::cout << "Setting root position to: " << _p << std::endl;
  _rootPosition = _p;
}

}
}

#endif /* _ROBOT_H_ */
