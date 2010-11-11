#include <stdio.h>
#include <stdlib.h>
#include <cstdarg>
#include <iostream>

#include <armadillo>

#include "robot/robot.h"
#include "robot/tree.h"

using namespace std;
using namespace arma;
using namespace edu_berkeley_cs184::robot;

template <class T>
static inline vector<T> makeVector(size_t count, ...) {
  vector<T> vec;
	va_list ap;
	va_start(ap, count);
	for (size_t i = 0; i < count; i++) 
		vec.push_back(va_arg(ap, T));
	va_end(ap);
	return vec;
}

static inline vec3 makeVec3(double x1, double x2, double x3) {
	vec3 v;
	v[0] = x1;
	v[1] = x2;
	v[2] = x3;
	return v;
}

int main(int argc, char **argv) {	

	TreeNode* root = new INode(
  	makeVector<LinkState*>(2, new LinkState(1.0, makeVec3(0, 0, 1), 0.0, makeVec3(-1, -1, 0)), new LinkState(1.0, makeVec3(0, 0, 1), 0.0, makeVec3(1, -1, 0))), 
		makeVector<TreeNode*>(2, new LNode(), new LNode()));

	LinkedTreeRobot robot(makeVec3(0, 0, 0), root);

	cout << "numJoints: " << robot.getNumJoints() << ", " <<
					"numEffectors: " << robot.getNumEffectors() << endl;

	mat J; 
	robot.computeJacobian(J);

	cout << J << endl;

  return 0;
}
