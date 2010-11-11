#include "robot.h"

#include <cassert>

using namespace arma;
using namespace std;

namespace edu_berkeley_cs184 {
namespace robot {

/**
 * Takes ownership of root
 */
LinkedTreeRobot::LinkedTreeRobot(const vec3& pos, TreeNode* root) : _rootPosition(pos), _root(root) {
  assert(_root != NULL && !root->isLeafNode());

  _root->assignNodeIndicies(0, 0); // assign indicies to each of the joints and effectors separately, index based from 0.

  // compute the number of joints & number of effectors
  _numJoints    = _root->numINodes();
  _numEffectors = _root->numLeafNodes();

  _root->gatherLeaves(_effectors);
  assert(_effectors.size() == _numEffectors);
}

LinkedTreeRobot::~LinkedTreeRobot() {
  delete _root;
}

mat& LinkedTreeRobot::computeJacobian(mat& m) const {
  m.zeros(_numEffectors * 3, _numJoints);

  // for each effector, traverse up its link hierarchy and compute
  // J(i,j) = (partial S_i / partial theta_j)

  for (vector<TreeNode*>::const_iterator it = _effectors.begin();
      it != _effectors.end();
      ++it) {

    TreeNode *effector = *it;
    size_t effectorId = effector->getIdentifier();

    vec3 effectorPos = effector->getGlobalPosition(_rootPosition);

    TreeNode *prevNode = effector;
    TreeNode *curNode = effector->getParent();
    while (curNode != NULL) {
      vec3 rotAxis = prevNode->getRotationAxis();
      size_t currentNodeId = curNode->getIdentifier();

      vec3 jacobianEntry = cross(rotAxis, effectorPos - curNode->getGlobalPosition(_rootPosition)); 
      m(3 * effectorId,     currentNodeId) = jacobianEntry[0];
      m(3 * effectorId + 1, currentNodeId) = jacobianEntry[1];
      m(3 * effectorId + 2, currentNodeId) = jacobianEntry[2];

      prevNode = curNode;
      curNode = curNode->getParent();
    }
  }

  return m;
}

}
}
