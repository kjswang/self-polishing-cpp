#ifndef _DISCON3D_H_
#define _DISCON3D_H_

/*
* Robot-obstacle distance calculation library 
* Author - 
* Date - 
*/

#include <vector>
#include <Eigen/Dense>
#include "structure.h"
#include "CapPos.h"
#include <cmath>
#include "math_pack.h"

using namespace std;
using namespace Eigen;

double distLinSeg2PC(MatrixXd P, MatrixXd Q, MatrixXd PC){
    MatrixXd PQ = P - Q;
    MatrixXd PC_Q = PC - Q;

    if (PQ.norm() == 0) {
        return PC_Q.colwise().norm().minCoeff();
    }
    MatrixXd T = PC_Q.transpose() * PQ / (PQ.transpose() * PQ);
    T = T.cwiseMax(0).cwiseMin(1);
    MatrixXd TQ = PQ.transpose() * T.transpose();
    MatrixXd dist = (PC_Q - TQ).colwise().norm();

    return dist.minCoeff();
}



double dist_arm_PC(MatrixXd theta, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, int check = 0){
    
    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;

    lineseg* pos = CapPos_real(base, DH, cap, tl);

    MatrixXd point1, point2;

    double dist, radius;
    vector<double> distList(sizeof(pos));
    for (int i = 0; i < nlink+2; i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC(point1, point2, PC);
        dist -= radius;

        distList[i] = dist;
    }

    d = min_element(distList.begin(),distList.end());
    return d;
}



double dist_arm_PC_WP(MatrixXd x0, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, int check = 0){
    MatrixXd theta, wp_pos;
    MatrixXd zero_tmp(1,1);
    zero_tmp(0,0) = 0;
    theta = x0.block(0,0,6,x0.cols());
    wp_pos = Vcat(zero_tmp,x0.row(6));
    wp_pos = Vcat(wp_pos,x0.row(7));

    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;

    lineseg* pos = CapPos_real(base, DH, cap, tl);

    MatrixXd point1, point2;
    double dist,radius;
    vector<double> distList(sizeof(pos));
    
    for (int i = 0; i < nlink+2; i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC(point1, point2, PC);
        dist -= radius;

        distList[i] = dist;
    }

    d = min_element(distList.begin(),distList.end());
    return d;
}

#endif