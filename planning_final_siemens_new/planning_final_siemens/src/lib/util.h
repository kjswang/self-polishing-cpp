#ifndef _CommonUtil_
#define _CommonUtil_

/*
* Common utility functinos for plannig cpp files
* Author - Weiye Zhao
* Date - Dec 29, 2019
*/

#include <Eigen/Dense>
#include "structure.h"
#include <cmath>

using namespace Eigen;
using namespace std;



void CapPos1(MatrixXd base, MatrixXd DH, capsule RoCap[], MatrixXd* M, lineseg* pos){
  int nlink = DH.rows();
  MatrixXd R, T, JTR;
  M[0].resize(4,4);
  M[0] << MatrixXd::Identity(3,3), base,
          0, 0, 0, 1;

  for (int i=0; i<nlink; ++i){
    R.resize(3,3);
    R << cos(DH(i,0)), -sin(DH(i,0))*cos(DH(i,3)), sin(DH(i,0))*sin(DH(i,3)),
          sin(DH(i,0)), cos(DH(i,0))*cos(DH(i,3)), -cos(DH(i,0))*sin(DH(i,3)),
          0, sin(DH(i,3)), cos(DH(i,3));

    T.resize(3,1);
    T << DH(i,2) * cos(DH(i,0)), 
         DH(i,2) * sin(DH(i,0)), 
         DH(i,1);

    JTR.resize(4,4);
    JTR << R, T,
           MatrixXd::Zero(1,3), 1;
       
    M[i+1].resize(4,4);
    M[i+1] = M[i]*JTR;

    // update the end-point position of capsule
    pos[i].p1 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(0) + M[i+1].block(0,3,3,1);
    pos[i].p2 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(1) + M[i+1].block(0,3,3,1);
  }
  }


#endif