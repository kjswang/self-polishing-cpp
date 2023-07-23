#ifndef _Structure_
#define _Structure_

/*
* CFS strucutre defintion 
* specify all the referenced strcture here
* Author -
* Date - 
*/

#include <Eigen/Dense>
#include <string>
#include <vector>


using namespace Eigen;
using namespace std;

struct capsule{
    MatrixXd p;
    double r;
};

struct tool{
	MatrixXd p1;
	MatrixXd p2;
	MatrixXd p3;
	double r1;
	double r2;
	double r3;
};

struct lineseg{
    MatrixXd p1;
    MatrixXd p2;
    double r;
};

struct inputfield{
	string field;
	vector<float> elements;
};

#endif