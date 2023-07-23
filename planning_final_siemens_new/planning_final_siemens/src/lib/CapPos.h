#ifndef _CapPos_
#define _CapPos_

/*
* Forward kinematics for capsules
* Author - Weiye Zhao
* Date - Nov. 5, 2019
*/

#include <Eigen/Dense>
#include "structure.h"
#include <cmath>
#include "load_param.h"
#include "global_var.h"

using namespace Eigen;
using namespace std;


void CapPos(MatrixXd base, MatrixXd DH, capsule RoCap[], MatrixXd* M, lineseg* pos){
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


lineseg* CapPos_new(MatrixXd base, MatrixXd DH, capsule RoCap[]){
    int nlink = DH.rows();
    lineseg* pos = new lineseg[nlink];
    MatrixXd M[nlink+1];
    MatrixXd R, T, JTR;

    M[0] = MatrixXd::Identity(4,4);
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
        pos[i].p1 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(0) + M[i+1].block(0,3,3,1) + base;
        pos[i].p2 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(1) + M[i+1].block(0,3,3,1) + base;
        pos[i].r = RoCap[i].r;
    }
    return pos;
}

lineseg* CapPos_real(MatrixXd base, MatrixXd DH, capsule RoCap[], tool tl){
    // MatrixXd Msix2tool;
    // loadjnt2tool(Msix2tool);
    int nlink = DH.rows();
    lineseg* pos = new lineseg[nlink+2];
    MatrixXd M[nlink+1];
    MatrixXd R, T, JTR;

    M[0] = MatrixXd::Identity(4,4);
    for (int i=0; i<nlink; ++i){
        if (i < nlink-1){
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
            pos[i].p1 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(0) + M[i+1].block(0,3,3,1) + base;
            pos[i].p2 = M[i+1].block(0,0,3,3) * RoCap[i].p.col(1) + M[i+1].block(0,3,3,1) + base;
            pos[i].r = RoCap[i].r;
        }
        if (i == nlink-1){
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
            // addtional matrix for the transformation from joint 6 to toolframe
            M[i+1] = M[i+1]*Msix2tool;
            // special design for fixture tool capsules 
            pos[i].p1 = M[i+1].block(0,0,3,3) * tl.p1.col(0) + M[i+1].block(0,3,3,1) + base;
            pos[i].p2 = M[i+1].block(0,0,3,3) * tl.p1.col(1) + M[i+1].block(0,3,3,1) + base;
            pos[i].r = tl.r1;
            // special design for fixture tool capsules 
            pos[i+1].p1 = M[i+1].block(0,0,3,3) * tl.p2.col(0) + M[i+1].block(0,3,3,1) + base;
            pos[i+1].p2 = M[i+1].block(0,0,3,3) * tl.p2.col(1) + M[i+1].block(0,3,3,1) + base;
            pos[i+1].r = tl.r2;
            // special design for fixture tool capsules 
            pos[i+2].p1 = M[i+1].block(0,0,3,3) * tl.p3.col(0) + M[i+1].block(0,3,3,1) + base;
            pos[i+2].p2 = M[i+1].block(0,0,3,3) * tl.p3.col(1) + M[i+1].block(0,3,3,1) + base;
            pos[i+2].r = tl.r3;
        }
    }
    return pos;
}


void ForKine_control(MatrixXd theta_ini, MatrixXd DH, MatrixXd base, MatrixXd& pos, MatrixXd& axisangle){
    /*
    * The Forward Kinematics to transform the original 
    * 6-dimensional angles data to end-effector pos 
    * and the quaternion q 
    */
    assert(theta_ini.rows()==6); // make sure it is a column vector
    // MatrixXd Msix2tool;
    // loadjnt2tool(Msix2tool);
    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta_ini;
    assert(DH.block(0,0,nlink,1) == theta_ini);
    MatrixXd M[nlink+1];
    MatrixXd R, T, JTR;
    M[0].resize(4,4);
    M[0] << MatrixXd::Identity(4,4);
    for (int i=0; i<nlink; ++i){
        // Rotation
        R.resize(3,3);
        R << cos(DH(i,0)), -sin(DH(i,0))*cos(DH(i,3)), sin(DH(i,0))*sin(DH(i,3)),
              sin(DH(i,0)), cos(DH(i,0))*cos(DH(i,3)), -cos(DH(i,0))*sin(DH(i,3)),
              0, sin(DH(i,3)), cos(DH(i,3));
        // Transformation
        T.resize(3,1);
        T << DH(i,2) * cos(DH(i,0)), 
             DH(i,2) * sin(DH(i,0)), 
             DH(i,1);

        JTR.resize(4,4);
        JTR << R, T,
               MatrixXd::Zero(1,3), 1;
           
        M[i+1].resize(4,4);
        M[i+1] = M[i]*JTR;

        // additional transformation matrix
        if (i == nlink-1){
             M[i+1]=M[i+1]*Msix2tool;
        }
    }
    // get the end-effector position 
    MatrixXd epos;
    epos = M[nlink].block(0,0,3,3)*MatrixXd::Zero(3,1) + M[nlink].block(0,3,3,1) + base;
    pos = epos;
    Matrix3d mat = M[nlink].block(0,0,3,3);
    // q = mat;
    AngleAxisd ang(mat);
    Vector3d agl;
    agl(0) = ang.axis()(0) * ang.angle();
    agl(1) = ang.axis()(1) * ang.angle();
    agl(2) = ang.axis()(2) * ang.angle();
    // Vector3d eagl = mat.eulerAngles(0,1,2);
    axisangle.resize(3,1);
    axisangle << agl(0),agl(1),agl(2);
}


#endif