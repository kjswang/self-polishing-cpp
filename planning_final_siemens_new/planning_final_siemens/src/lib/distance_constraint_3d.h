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

double distLinSeg(MatrixXd s1, MatrixXd e1, MatrixXd s2, MatrixXd e2){
  MatrixXd d1, d2, d12;
  d1 = e1 - s1;
  d2 = e2 - s2;
  d12 = s2 - s1;
  
  double D1 = pow_sum(d1, 2);
  double D2 = pow_sum(d2, 2);
  double S1 = dot_sum(d1, d12);
  double S2 = dot_sum(d2, d12);
  double R = dot_sum(d1, d2);
  double den = D1*D2 - pow(R,2);


  //if one of the segments is a point
  double u,t,uf;
  if (D1 == 0 || D2 == 0){
    if (D1 != 0){
      u = 0;
      t = S1 / D1;
      t = fixbound(t);
    }
    else if (D2 != 0){ 
      t = 0;
      u = -S2 / D2;
      u = fixbound(u);
    }
    else{ // both line segment is point 
      t = 0;
      u = 0;
    }
  }

  else if (den == 0){ // line is parallel 
    t = 0;
    u = -S2/D2;
    uf = fixbound(u);
    if (uf != u){
      t = (uf*R+S1)/D1;
      t = fixbound(t);
      u = uf;
    }
  }

  else{ // general case
    t = (S1*D2-S2*R)/den;
    t = fixbound(t);
    u = (t*R-S2)/D2;
    uf = fixbound(u);
    if (uf != u){
      t = (uf*R+S1)/D1;
      t = fixbound(t);
      u = uf;
    }
  }

    // compute the distance, given t and u 
    MatrixXd distm;
    distm = d1*t-d2*u-d12;
    double dist = distm.norm();
    // point
    MatrixXd point(3,2);
    point.block(0,0,3,1) = s1 + d1*t;
    point.block(0,1,3,1) = s2 + d2*u;
    return dist;

}


double dist_arm_3D_Heu(MatrixXd& theta, MatrixXd& DH, MatrixXd& base, MatrixXd& obs, capsule cap[]){
  int nstate = theta.rows();
  for (int i=0; i<nstate; ++i){
    DH(i,0) = theta(i,0);
  }
  double d = std::numeric_limits<double>::infinity(); //set d as positive infinity

  // forward kinematics for position calculation under theta
  lineseg pos[DH.rows()];
  MatrixXd M[DH.rows()+1];
  CapPos(base, DH, cap, M, pos);
  // cout << base << endl << DH << endl;
  
  // calculate the closest distance 
  double dis;
  int link_id;
  for (int i=0; i<nstate; ++i){
    dis = distLinSeg(pos[i].p1, pos[i].p2, obs.block(0,0,3,1), obs.block(0,1,3,1)); // this is 3d scenario. where obstacle position in 3d form
    if (dis < d){
      d = dis;
      link_id = i+1;
    }
  }
  return d;
}


double distLinSeg2PC(MatrixXd P, MatrixXd Q, MatrixXd PC){
    MatrixXd PQ = P - Q;
    MatrixXd PC_Q = PC - Q;

    if (PQ.norm() == 0) {
        return PC_Q.colwise().norm().minCoeff();
    }

    double t = (PC_Q.transpose()*PQ).norm() / (PQ.transpose()*PQ).norm();
    t = max(0.0, min(1.0,t)); 

    MatrixXd TQ = t * PQ.transpose();
    MatrixXd dist = (PC_Q - TQ).colwise().norm();

    // MatrixXd T = PC_Q.transpose() * PQ / (PQ.transpose() * PQ);
    // T = T.cwiseMax(0).cwiseMin(1);
    // MatrixXd TQ = PQ.transpose() * T.transpose();
    // MatrixXd dist = (PC_Q - TQ).colwise().norm();

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

    d = *std::min_element(distList.begin(),distList.end());
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

    d = *std::min_element(distList.begin(),distList.end());
    return d;
}

#endif