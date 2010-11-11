#include "tree.h"

#include <stdexcept>
#include <cassert>
#include <stack>

using namespace std;
using namespace arma;

namespace edu_berkeley_cs184 {
namespace robot {

TreeNode::TreeNode() { _parent = NULL; }
TreeNode::~TreeNode() {}

vec3 TreeNode::getGlobalPosition(const vec3& rootPoint) {
  stack<TreeNode*> pathToRoot;
  stack<size_t> indexes;
  TreeNode* cur = getParent();
  TreeNode* prev = this;
  while (cur != NULL) {
    pathToRoot.push(cur);
    indexes.push(prev->getIndex());
    prev = cur;
    cur = cur->getParent();
  }

  assert(pathToRoot.size() == indexes.size());

  vec3 curPos = rootPoint;
  while (!pathToRoot.empty()) {
    TreeNode *n = pathToRoot.top();
    size_t idx = indexes.top();
    pathToRoot.pop(); indexes.pop();
    vector<LinkState*>::const_iterator state_iter = n->getLinkStates() + idx;
    curPos = (*state_iter)->getEndpoint(curPos);
  }
  return curPos;
}


/**
 * Takes ownership of objs in states and objs in children
 */
INode::INode(const vector<LinkState*>& states, 
             const vector<TreeNode*>& children) : _states(states), _kids(children) {
  assert(states.size() == children.size());
  assert(states.size() > 0); // INodes are required to have >= 1 child
  for (size_t i = 0; i < _kids.size(); i++) {
    _kids[i]->setParent(this);
    _kids[i]->setIndex(i);
  }
}

/**
 * Recursively deletes link states first, and then children nodes
 */
INode::~INode() {
  for (vector<LinkState*>::iterator it = _states.begin();
      it != _states.end();
      ++it)
    delete *it;
  for (vector<TreeNode*>::iterator it = _kids.begin();
      it != _kids.end();
      ++it)
    delete *it;
}

bool INode::isLeafNode() const { return false; }

vector<TreeNode*>::const_iterator INode::getChildren() const { return _kids.begin(); }
vector<LinkState*>::const_iterator INode::getLinkStates() const { return _states.begin(); }

size_t INode::numLeafNodes() const {
  size_t sumSoFar = 0;
  for (vector<TreeNode*>::const_iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    sumSoFar += (*it)->numLeafNodes();
  return sumSoFar;
}

size_t INode::numINodes() const {
  size_t sumSoFar = 1;
  for (vector<TreeNode*>::const_iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    sumSoFar += (*it)->numINodes();
  return sumSoFar;
}

pair<size_t, size_t> INode::assignNodeIndicies(const size_t iNodeOffset, const size_t lNodeOffset) {
  // must do this in prefix order.
  idx = iNodeOffset;
  size_t numINodes = 1;
  size_t numLNodes = 0;
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it) {
    pair<size_t, size_t> thisNode = (*it)->assignNodeIndicies(iNodeOffset + numINodes, 
                                                              lNodeOffset + numLNodes);
    numINodes += thisNode.first;
    numLNodes += thisNode.second;
  }
  return pair<size_t, size_t>(numINodes, numLNodes);
}

vector<TreeNode*>& INode::gatherLeaves(vector<TreeNode*>& buffer) {
  for (vector<TreeNode*>::iterator it = _kids.begin(); 
      it != _kids.end();
      ++it)
    (*it)->gatherLeaves(buffer);
  return buffer;
}


bool LNode::isLeafNode() const { return true; }

vector<TreeNode*>::const_iterator LNode::getChildren() const { 
  throw runtime_error("getChildren on leaf node"); 
}

vector<LinkState*>::const_iterator LNode::getLinkStates() const { 
  throw runtime_error("getLinkStates on leaf node"); 
}

size_t LNode::numLeafNodes() const { return 1; }

size_t LNode::numINodes() const { return 0; }

pair<size_t, size_t> LNode::assignNodeIndicies(const size_t iNodeOffset, const size_t lNodeOffset) {
  idx = lNodeOffset;
  return pair<size_t, size_t>(0, 1);
}

vector<TreeNode*>& LNode::gatherLeaves(vector<TreeNode*>& buffer) {
  buffer.push_back(this);
  return buffer;
}

}
}
