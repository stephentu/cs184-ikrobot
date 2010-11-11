#ifndef _ROBOT_H_
#define _ROBOT_H_

#include <armadillo>

#include "tree.h"

namespace edu_berkeley_cs184 {
namespace robot {

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

  arma::mat& computeJacobian(arma::mat&) const;
  arma::vec computeDeltaThetas(const arma::vec&) const;
  void updateThetas(const arma::vec&);

  /** Main IK method- given an input of desired positions, update the robot to
   * try to match this input */
  void solveIK(const arma::vec&);

  /** draw the robot in openGL */
  void renderRobot() const;

private:
  const arma::vec3 _rootPosition;
  TreeNode* _root;

  // computed from _root
  size_t _numJoints;
  size_t _numEffectors;

  // effectors, from traversing the root
  std::vector<TreeNode*> _effectors;
};

inline size_t LinkedTreeRobot::getNumJoints() const {
  return _numJoints;
}

inline size_t LinkedTreeRobot::getNumEffectors() const {
  return _numEffectors;
}

}
}

#endif /* _ROBOT_H_ */
