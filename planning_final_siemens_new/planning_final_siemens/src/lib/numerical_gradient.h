#ifndef _NGrad_
#define _NGrad_

/*
* Numerical gradient library 
* Use 2nd order central difference 
* Author - 
* Date - 
*/

#include <Eigen/Dense>
#include "structure.h"
#include "distance_constraint_3d.h"
#include <cmath>
#include "safe_polishing.h"

using namespace Eigen;
using namespace std;


MatrixXd central_diff(MatrixXd theta, MatrixXd DH, MatrixXd base, MatrixXd obs, capsule cap[]){
    double dist = dist_arm_3D_Heu(theta, DH, base, obs, cap);
    MatrixXd grad;

    grad.resize(1, theta.rows());
    double eps = 1e-5;

    MatrixXd theta_tmp;
    theta_tmp = theta;
    for (int i=0; i<theta.rows(); ++i){
        theta_tmp(i,0) = theta(i,0) + eps/2;
        double dist_h = dist_arm_3D_Heu(theta_tmp, DH, base, obs, cap);
        theta_tmp(i,0) = theta(i,0) - eps/2;
        double dist_l = dist_arm_3D_Heu(theta_tmp, DH, base, obs, cap);
        grad(0,i) = (dist_h - dist_l) / eps;
    }
    // cout << grad;
    return grad;
}


MatrixXd central_diff_polish(MatrixXd theta, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC){
    int dim = theta.rows();
    double y = dist_arm_PC(theta, DH, base, cap, tl, PC);
    MatrixXd grad(1, dim);

    double eps = 1e-5;
    MatrixXd x = theta;
    MatrixXd xp = theta;
    assert(theta.cols() == 1);

    double yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = dist_arm_PC(theta, DH, base, cap, tl, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = dist_arm_PC(theta, DH, base, cap, tl, PC);
        grad(0,i) = (yhi - ylo)/eps;
        xp(i,0) = x(i,0);
    }
    return grad; 
}


MatrixXd central_diff_polish_collaboration(MatrixXd x0, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC){
    int dim = x0.rows();
    double y = dist_arm_PC_WP(x0, DH, base, cap, tl, PC);
    MatrixXd grad(1, dim);

    double eps = 1e-5;
    MatrixXd x = x0;
    MatrixXd xp = x0;
    assert(x0.cols() == 1);

    double yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = dist_arm_PC_WP(x0, DH, base, cap, tl, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = dist_arm_PC_WP(x0, DH, base, cap, tl, PC);
        grad(0,i) = (yhi - ylo)/eps;
        xp(i,0) = x(i,0);
    }
    return grad; 
}

MatrixXd Diff_Jac_num_grad(MatrixXd x0, MatrixXd c1, MatrixXd base, capsule cap[], tool tl, MatrixXd PC){
    int dim = x0.rows();
    MatrixXd y = next_point_WP(x0, c1, PC);
    MatrixXd grad(1, dim);

    double eps = 1e-5;
    MatrixXd x = x0;
    MatrixXd xp = x0;
    assert(x0.cols() == 1);

    MatrixXd yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = next_point_WP(x0, c1, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = next_point_WP(x0, c1, PC);
        grad(0,i) = (yhi - ylo).value()/eps;
        xp(i,0) = x(i,0);
    }
    return grad; 
}

#endif