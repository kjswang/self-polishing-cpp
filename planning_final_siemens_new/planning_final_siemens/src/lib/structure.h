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
#include <armadillo>

using namespace Eigen;
using namespace std;
using namespace arma;

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

struct capsule_arma{
    mat p;
    double r;
};

struct tool_arma{
	mat p1;
	mat p2;
	mat p3;
	double r1;
	double r2;
	double r3;
};

struct lineseg_arma{
    mat p1;
    mat p2;
    double r;
};


struct inputfield{
	string field;
	vector<float> elements;
};


void loadDHbase_arma(mat& DH, mat& base){
    DH << 0 << 0.281 << 0.145 << -1.570796 << endr
           << -1.570796 << 0 << 0.87 << 0 << endr
           << 0 << 0 << 0.21 << -1.570796 << endr
           << 0 << 1.025 << 0 << 1.570796 << endr
           << 0 << 0 << 0 << -1.570796 << endr
           << 0 << 0.175 << 0 << 0 << endr;
    base << 0 << endr
            << 0 << endr
            << 0.259 << endr;
}

void loadCapsules(capsule_arma* cap) {
    // Load information for each capsule
    // 0
    cap[0].p.resize(3, 2);
    cap[0].p << -0.145 << -0.145 << endr
              <<  0.105 <<  0.105 << endr
              <<  0 << 0 << endr;
    cap[0].r = 0.385;
    // 1
    cap[1].p.resize(3, 2);
    cap[1].p << -0.87 << 0 << endr
              <<  0 <<  0 << endr
              <<  -0.1945 <<  -0.1945 << endr;
    cap[1].r = 0.195;
    // 2
    cap[2].p.resize(3, 2);
    cap[2].p << -0.02 <<  -0.09 << endr
             <<   0.073 << 0.073 << endr
             <<   0.115 <<  0.115 << endr;
    cap[2].r = 0.33;
    // 3
    cap[3].p.resize(3, 2);
    cap[3].p << 0 <<  0 << endr
             <<   -0.65 <<  0 << endr
             <<   -0.0235 <<  -0.0235 << endr;
    cap[3].r = 0.115;
    // 4
    cap[4].p.resize(3, 2);
    cap[4].p << 0 << 0 << endr
             <<   0.0145 << 0.0145 << endr
            << 0.025 << 0.025 << endr;
    cap[4].r = 0.15;
}


void loadTool_arma(tool_arma& tl) {
    tl.p1.resize(3, 2);
    tl.p1 << -0.142391115489123 << -0.00817305105262633 << endr
            << 0.0844531341176073 <<  0.0515485707207726 << endr
             << -0.448552826957168 <<  -0.216598563316495 << endr;

    tl.p2.resize(3, 2);
    tl.p2 << -0.00407589262314147 << -0.00777057759414465 << endr
            << 0.0819148210462628 << 0.0218344612710158 << endr
              << -0.209517946591771 << -0.215902972600098 << endr;

    mat pleft, pright;
    pleft.resize(3, 1);
    pright.resize(3, 1);
    pleft << 0.00252978728478826 << endr 
            << 6.28378116607958e-10 << endr 
            << -0.170767309373314 << endr;
    pright << 0.000390496336481267 << endr 
            << 1.00300828261106e-10 << endr 
            << -0.0344384157898974 << endr;
    mat diff = pright - pleft;
    mat prnew = pleft + 0.3 * diff;

    tl.p3 = join_horiz(pleft, prnew);
    tl.r1 = 0.065;
    tl.r2 = 0.05;
    tl.r3 = 0.03;
}




#endif