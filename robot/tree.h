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

  /** Assign indicies in prefix order.
   * 1st param- node offset
   * 2nd param- joint offset
   * 3rd param- leaf node offset
   * 4th param- num nodes set
   * 5th param- num joint nodes set
   * 6th param- num leaf nodes set
   */ 
  virtual void assignNodeIndicies(const size_t, 
                                  const size_t,
                                  const size_t,
                                  size_t&,
                                  size_t&,
                                  size_t&) = 0;

  /** Gather all the leaves of this tree into buffer, using push_back to
   * append */
  virtual std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&) = 0;

  /** Gather all inner nodes of this tree (including this node) into buffer,
   * using push_back to insert */
  virtual std::vector<TreeNode*>& gatherInnerNodes(std::vector<TreeNode*>&) = 0;

  /** Gather all nodes of this tree (including this node) into buffer,
   * using push_back to insert */
  virtual std::vector<TreeNode*>& gatherNodes(std::vector<TreeNode*>&) = 0;

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

  /** An identifier for every node. */
  inline size_t getIdentifier() const;

  /** An identifier for leaf nodes only. is an error for non leaf nodes */
  virtual size_t getLeafIdentifier() const = 0;

  /** Update all angle by being given deltas. rotation axes given as a
   * reference */
  virtual void updateThetas(const arma::vec&, const arma::vec&) = 0;

  virtual void renderTree(Context&) const = 0;

  /** builds a context for THIS node */
  Context& getContextForNode(Context&);

  inline bool isFixed() const;

  inline void setFixed(const bool);

  inline arma::vec3 getFixedPosition() const;

  inline void setFixedPosition(const arma::vec3&);

protected:
  TreeNode* _parent;
  size_t idx;
  size_t identity;
  bool _fixed; // is there a positional constraint
  arma::vec3 _position; // what is the fixed position constraint?
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

inline size_t TreeNode::getIdentifier() const { return identity; }

inline bool TreeNode::isFixed() const { return _fixed; }

inline void TreeNode::setFixed(const bool _f) { _fixed = _f; }

inline arma::vec3 TreeNode::getFixedPosition() const { return _position; }

inline void TreeNode::setFixedPosition(const arma::vec3& _p) { _position = _p; }

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
  void assignNodeIndicies(const size_t, 
                          const size_t,
                          const size_t,
                          size_t&,
                          size_t&,
                          size_t&);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&);
  std::vector<TreeNode*>& gatherInnerNodes(std::vector<TreeNode*>&);
  std::vector<TreeNode*>& gatherNodes(std::vector<TreeNode*>&);
  void updateThetas(const arma::vec&, const arma::vec&);
  void renderTree(Context&) const;
  size_t getLeafIdentifier() const;
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
  void assignNodeIndicies(const size_t, 
                          const size_t,
                          const size_t,
                          size_t&,
                          size_t&,
                          size_t&);
  std::vector<TreeNode*>& gatherLeaves(std::vector<TreeNode*>&);
  std::vector<TreeNode*>& gatherInnerNodes(std::vector<TreeNode*>&);
  std::vector<TreeNode*>& gatherNodes(std::vector<TreeNode*>&);
  void updateThetas(const arma::vec&, const arma::vec&);
  void renderTree(Context&) const;
  size_t getLeafIdentifier() const;
private:
  size_t leafIdentity;
};

}
}

#endif /* _TREE_H_ */
