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
  _numJoints    = _root->numEdges();
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
      size_t jointId = prevNode->getEdgeIdentifier();

      vec3 jacobianEntry = cross(rotAxis, effectorPos - curNode->getGlobalPosition(_rootPosition)); 
      m(3 * effectorId,     jointId) = jacobianEntry[0];
      m(3 * effectorId + 1, jointId) = jacobianEntry[1];
      m(3 * effectorId + 2, jointId) = jacobianEntry[2];

      prevNode = curNode;
      curNode = curNode->getParent();
    }
  }

  return m;
}

}
}
