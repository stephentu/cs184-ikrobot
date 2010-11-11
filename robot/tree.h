#ifndef _TREE_H_
#define _TREE_H_

#include <vector>
#include <utility>
#include <armadillo>
#include <cassert>

#include "../util/util.h"

namespace edu_berkeley_cs184 {
namespace robot {

/**
 * Mutable link state
 */
class LinkState {
public:
  LinkState(const double _l,
            const arma::vec3& _axis,
            const double _t0,
            const arma::vec3& _baseline) :
    length(_l), axis(_axis), angle(_t0), baselineDirection(_baseline) {

    assert(_l > 0.0);

    // normalize axis and baseline just to be sure
    edu_berkeley_cs184::util::normalize_vec3(axis); 
    edu_berkeley_cs184::util::normalize_vec3(baselineDirection);

    // now assert that the dot product of axis and baselineDirection is zero
    assert( edu_berkeley_cs184::util::double_equals(arma::dot(axis, baselineDirection), 0) );
  }

  /** Given an input coordinate, compute the endpoint of this link given the
   * jointLoc (the coordinate of the joint) */
  inline arma::vec3 getEndpoint(const arma::vec3& jointLoc) const {
    // step 1: rotate baselineDirection by angle degrees along axis vector (in
    // right hand rule sense)
    arma::vec3 rotated = edu_berkeley_cs184::util::rotate_expmap(baselineDirection, axis, angle);

    // step 2: endpoint = jointLoc + length * rotated
    return jointLoc + length * rotated;
  }

  const double length; /* not sure what unit this is */
  arma::vec3 axis; /* NORMALIZED axis of rotation */ 
  double angle; /* radians */
  arma::vec3 baselineDirection; /* When angle = 0, what NORMALIZED dir is this link pointing in? */
};

class TreeNode {
public:
  TreeNode();
  virtual ~TreeNode();

  virtual bool isLeafNode() const = 0;

  virtual std::vector<TreeNode*>::const_iterator getChildren() const = 0;
  virtual std::vector<LinkState*>::const_iterator getLinkStates() const = 0;

  virtual size_t numLeafNodes() const = 0;
  virtual size_t numINodes() const = 0;

  /** assign node indicies in prefix order for inodes and leaf nodes
   * separately, using the inputs as the beginning offset. returns
   * the number of intermediate and leaf nodes assigned respectively */ 
  virtual std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t) = 0;

  /** Gather all the leaves of this tree into buffer, using push_back to
   * append */
  virtual std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>& ) = 0;

  /** Returns the global position of THIS node, given the input position for
   * the ROOT node. If this node IS the root node, then the same position will
   * be returned. */
  arma::vec3 getGlobalPosition(const arma::vec3&);

  inline bool isRootNode() const;

  inline TreeNode* getParent() const;

  /** TODO: take out of public access */
  inline void setParent(TreeNode*);

  /** TODO: take out of public access */
  inline void setIndex(const size_t);

  /** The index (0-based) to look in the parent's child array for me */
  inline size_t getIndex() const;

  /** An identifier assigned to me based on my position in the entire tree.
   * DO NOT GET THIS CONFUSED WITH getIndex(). assignNodeIndicies MUST be
   * called BEFORE this method is called for this ID to be meaningful */
  inline size_t getIdentifier() const;

  /** What axis of rotation does THIS node rotate about (requires looking into
   * parent). is ill-defined for the root (since the root is fixed!) */
  inline arma::vec3 getRotationAxis() const; 

protected:
  TreeNode* _parent;
  size_t idx;
  size_t id;
};

inline bool TreeNode::isRootNode() const { return _parent == NULL; }

inline TreeNode* TreeNode::getParent() const { return _parent; }

inline void TreeNode::setParent(TreeNode *p) { _parent = p; }

inline void TreeNode::setIndex(const size_t _i) { idx = _i; }

inline size_t TreeNode::getIndex() const {
  assert(!isRootNode());
  return idx;
}

inline size_t TreeNode::getIdentifier() const { return id; }

inline arma::vec3 TreeNode::getRotationAxis() const {
  assert(!isRootNode());
  return (*(getParent()->getLinkStates() + getIndex()))->axis;
}


class INode : public TreeNode {
public:
  INode(const std::vector<LinkState*>&, const std::vector<TreeNode*>&);
  virtual ~INode();
  bool isLeafNode() const;
  std::vector<TreeNode*>::const_iterator getChildren() const;
  std::vector<LinkState*>::const_iterator getLinkStates() const;
  size_t numLeafNodes() const;
  size_t numINodes() const;
  std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>& );
private:
  std::vector<LinkState*> _states;
  std::vector<TreeNode*> _kids;
};

class LNode : public TreeNode {
public:
  bool isLeafNode() const;
  std::vector<TreeNode*>::const_iterator getChildren() const;
  std::vector<LinkState*>::const_iterator getLinkStates() const;
  size_t numLeafNodes() const;
  size_t numINodes() const;
  std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>& );
};

}
}

#endif /* _TREE_H_ */