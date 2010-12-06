#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include <stack>
#include <armadillo>

namespace edu_berkeley_cs184 {
namespace robot {

/** maintains context for rotations */
class Context {
public:
  // position (0, 0, 0), transform eye(4)
  Context();

  Context(const arma::vec3&);

  Context(const arma::vec3&,
          const arma::mat44&);

  void pushContext(const arma::vec3&,
                   const arma::vec3&,
                   const double);

  void pushContext(const arma::vec3&, 
                   const double,
                   const double,
                   const double);

  void pushContext(const arma::vec3&, const arma::mat44&);

	void pushContext(const arma::vec3&);

  void popContext();

  /** returns the current GLOBAL origin */
  inline arma::vec3 getCurrentOrigin() const;

  /** returns the given vector as a GLOBAL vector, assuming it is purely
   * directional (homogenous coordinate = 0) */
  inline arma::vec3 getVectorInContext(const arma::vec3&) const; 

private:
  std::stack<arma::vec3> _positions; // maintains the stack of positions
  std::stack<arma::mat44> _transforms; // maintains the stack of transformation matrices

};

inline arma::vec3 Context::getCurrentOrigin() const { return _positions.top(); }

inline arma::vec3 Context::getVectorInContext(const arma::vec3& v) const {
  return _transforms.top().submat(arma::span(0, 2), arma::span(0, 2)) * v;
}

}
}

#endif /* _CONTEXT_H */
