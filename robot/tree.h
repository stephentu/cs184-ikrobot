#ifndef _TREE_H_
#define _TREE_H_

#include <vector>
#include <utility>
#include <armadillo>
#include <cassert>

#include "../util/util.h"
#include "context.h"
#include "link.h"

namespace edu_berkeley_cs184 {
namespace robot {

class TreeNode {
public:
  TreeNode();
  virtual ~TreeNode();

  virtual bool isLeafNode() const = 0;

  virtual std::vector<TreeNode*>::const_iterator getChildren() const = 0;
  virtual std::vector<LinkState*>::const_iterator getLinkStates() const = 0;

  /** Number of leaf nodes at the bottom of tree */
  virtual size_t numLeafNodes() const = 0;

  /** Number of edges total in the tree (paths) */
  virtual size_t numEdges() const = 0;

  /** Number of DOF in the tree (by adding up all the link states) */
  virtual size_t numDOF() const = 0;

  /** assign node indicies in prefix order for inodes and leaf nodes
   * separately, using the inputs as the beginning offset. returns
   * the number of intermediate and leaf nodes assigned respectively */ 
  virtual std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t) = 0;

  /** Gather all the leaves of this tree into buffer, using push_back to
   * append */
  virtual std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&) = 0;

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

  /** My link state (requires lookup into parent, illdefined for root) */
  inline LinkState* getLinkState() const;

  /** An identifier for LEAF NODES only. Errors out for intermediate nodes */
  virtual size_t getIdentifier() const = 0;

  /** Update all angle by being given deltas. rotation axes given as a
   * reference */
  virtual void updateThetas(const arma::vec&, const arma::vec&) = 0;

  virtual void renderTree(Context&) const = 0;

  /** builds a context for THIS node */
  Context& getContextForNode(Context&);

protected:
  TreeNode* _parent;
  size_t idx;
};

inline bool TreeNode::isRootNode() const { return _parent == NULL; }

inline TreeNode* TreeNode::getParent() const { return _parent; }

inline void TreeNode::setParent(TreeNode *p) { _parent = p; }

inline void TreeNode::setIndex(const size_t _i) { idx = _i; }

inline size_t TreeNode::getIndex() const {
  assert(!isRootNode());
  return idx;
}

inline LinkState* TreeNode::getLinkState() const {
  assert(!isRootNode());
  return *(getParent()->getLinkStates() + getIndex());
}

class INode : public TreeNode {
public:
  INode(const std::vector<LinkState*>&, const std::vector<TreeNode*>&);
  virtual ~INode();
  bool isLeafNode() const;
  std::vector<TreeNode*>::const_iterator getChildren() const;
  std::vector<LinkState*>::const_iterator getLinkStates() const;
  size_t numLeafNodes() const;
  size_t numEdges() const;
  size_t numDOF() const;
  std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&);
  size_t getIdentifier() const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void renderTree(Context&) const;
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
  size_t numEdges() const;
  size_t numDOF() const;
  std::pair<size_t, size_t> assignNodeIndicies(const size_t, const size_t);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&);
  size_t getIdentifier() const;
  void updateThetas(const arma::vec&, const arma::vec&);
  void renderTree(Context&) const;
private:
  size_t id;
};

}
}

#endif /* _TREE_H_ */
