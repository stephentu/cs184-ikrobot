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

  arma::mat& computeJacobian(arma::mat&) const;
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
