#ifndef _ROBOT_H_
#define _ROBOT_H_

/*
* Robot property class defintion 
* robot property specification
* Author - 
* Date -
*/

#include <math.h>
#include <vector>
#include "math_pack.h"
#include <Eigen/Dense>
#include "CapPos.h"
#include "structure.h"
#include "global_var.h"
#include "load_param.h"


using namespace std;
using namespace Eigen;


class Robot {

    public:
    int nlink;
    int umax;
    double margin;
    double pi;
    double delta_t;
    string name;
    MatrixXd thetamax;
    MatrixXd lb; // upper joint angle bound for each joint 
    MatrixXd ub; // lower joint angle bound for each joint 
    MatrixXd thetadotmax;
    MatrixXd l;
    MatrixXd DH;
    MatrixXd base;
    MatrixXd A;
    MatrixXd B;
    MatrixXd Ac;
    MatrixXd Bc;
    MatrixXd goal; // the Cartesian position that end effector of robot should follow
    capsule cap[5];
    tool tl;
    lineseg pos[6];
    MatrixXd M[7];



  Robot(string robot_name){
    int robot_num = parse_name(robot_name);
    switch (robot_num){
        case 1: {
            M16iBproperty();
            break;
        }
        case 2: {
            GP50property();
            break;
        }
    }
    // initialize the kinematic matrix
    CapPos(base, DH, cap, M, pos);
  }


    int parse_name(string name){
    string gp50, M16iB;
    gp50 = "GP50";
    M16iB = "M16iB";
    if((name.compare(gp50)) == 0){
        cout << "robot is identified as " << gp50 << endl; 
        cout << "------------ start initialize GP50 property ------------" << endl;
        return 2;
    }

    if((name.compare(M16iB)) == 0){
        cout << "robot is identified as " << M16iB << endl; 
        cout << "------------ start initialize M16iB property ------------" << endl;
        return 1;
    }

    else{
        cout << "no match robot property found" << endl;
        abort();
    } 
    return 0;
  }


  void M16iBproperty(){
    name = "M16iB";
    pi = M_PI;
    nlink = 6;
    delta_t = 0.5;
    thetamax.resize(6,2);
    thetamax << -pi, pi,
                0, pi,
                -pi, pi,
                -pi, pi,
                -pi/2, pi/2,
                -pi, pi;

    thetadotmax.resize(6,1);
    thetadotmax << 1,
                   1,
                   1,
                   1,
                   0.5,
                   0.5;

    DH.resize(6,4);
    DH << 0.5, 0.65, 0.15, 1.5708,
          1.5708, 0, 0.77, 0,
          0, 0, 0.1, 1.5708,
          0, 0.74, 0, -1.5708,
          -pi/2, 0, 0, 1.5708,
          pi, 0.1, 0, 0;

    base.resize(3,1);
    base << 0,  
            0,
            0;


    // initialize the capsule of robot
    // 0
    cap[0].p.resize(3,2);
    cap[0].p << 0, 0,
                0, 0,
                -0.1, 0.1;      
    cap[0].r = 0.15;
    // 1
    cap[1].p.resize(3,2);
    cap[1].p << -0.75, 0,
                0, 0,
                -0.15, -0.15;    
    cap[1].r = 0.13;
    // 2
    cap[2].p.resize(3,2);
    cap[2].p << -0.03, -0.03,
                0, 0,
                0.05, 0.05;   
    cap[2].r = 0.22;
    // 3
    cap[3].p.resize(3,2);
    cap[3].p << 0, 0,
                0, 0.55,
                0, 0;
    cap[3].r = 0.11;
    // 4
    cap[4].p.resize(3,2);
    cap[4].p << 0, 0,
                0, 0,
                -0.05, 0.110;
    cap[4].r = 0.07;
    // 5
    cap[5].p.resize(3,2);
    cap[5].p << -0.11, -0.11,
                0, 0,
                0.09, 0.09;
    cap[5].r = 0.07;

    A.resize(12,12);
    A << MatrixXd::Identity(6,6), MatrixXd::Identity(6,6) * delta_t,
         MatrixXd::Zero(6,6),    MatrixXd::Identity(6,6);

    B.resize(12,6);
    B << MatrixXd::Identity(6,6) * 0.5 * pow(delta_t, 2.0),
         MatrixXd::Identity(6,6) * delta_t;

    Ac.resize(6,6);
    Ac << MatrixXd::Identity(3,3), MatrixXd::Identity(3,3) * delta_t,
          MatrixXd::Zero(3,3), MatrixXd::Identity(3,3);

    Bc.resize(6,3);
    Bc << MatrixXd::Identity(3,3) * 0.5 * pow(delta_t, 2.0),
          MatrixXd::Identity(3,3) * delta_t;
    }


    void GP50property(){
    name = "gp50";
    pi = M_PI;
    nlink = 6;
    delta_t = 0.5;
    thetamax.resize(6,2);
    thetamax << -pi, pi,
                -pi, pi,
                -pi, pi,
                -pi, pi,
                -pi, pi,
                -pi, pi;

    /* GP50 joints bounds */
    lb.resize(6,1);
    lb << -pi, 
          -pi,
          -pi, 
          -pi,
          -pi,
          -pi;
    ub.resize(6,1);
    ub << pi, 
          pi,
          pi, 
          pi,
          pi,
          pi;

    /* GP50 DH parameters and the base xyz position  */
    loadDHbase(DH, base);

    // initialize the capsule parameter of GP50
    // 0
    cap[0].p.resize(3,2);
    cap[0].p << -0.145, -0.145,
                0.105, 0.105,
                0, 0;    
    cap[0].r = 0.385;
    // 1
    cap[1].p.resize(3,2);
    cap[1].p << -0.87, 0,
                0, 0,
                -0.1945, -0.1945;    
    cap[1].r = 0.195;
    // 2
    cap[2].p.resize(3,2);
    cap[2].p << -0.02, -0.09,
                0.073, 0.073,
                0.115, 0.115;   
    cap[2].r = 0.33;
    // 3
    cap[3].p.resize(3,2);
    cap[3].p << 0, 0,
                -0.65, 0,
                -0.0235, -0.0235;
    cap[3].r = 0.115;
    // 4
    cap[4].p.resize(3,2);
    cap[4].p << 0, 0,
                0.0145, 0.0145,
                0.025, 0.025;
    cap[4].r = 0.15;

    // tool capusles are defined using tool strucutre 
    // tool is approximated using three capsules
    tl.p1.resize(3,2);
    tl.p1 << -0.142391115489123,-0.00817305105262633,
             0.0844531341176073,0.0515485707207726,
             -0.448552826957168,-0.216598563316495;
    tl.p2.resize(3,2);
    tl.p2 << -0.00407589262314147,-0.00777057759414465,
             0.0819148210462628,0.0218344612710158,
             -0.209517946591771,-0.215902972600098;

    // the last capsule will contact with weld-point 
    // only consider part of it;
    // the last capsule is shorter than real polishing tool length
    // this is compromised for satisfying conservative approximation constriants of workpiece 
    MatrixXd pleft, pright;
    double ratio;
    if (planner == 0){
        ratio = 1; // the ratio of last capsule in consideration when measurement planning 
    }
    if (planner == 1){
        // note tip will in contact with weld, only consider collision of part of it
        ratio = 0.3; // the ratio of last capsule in consideration when polishing planning
    }
    
    pleft.resize(3,1);
    pright.resize(3,1);
    pleft << 0.00252978728478826,6.28378116607958e-10,-0.170767309373314;
    pright << 0.000390496336481267,1.00300828261106e-10,-0.0344384157898974;
    MatrixXd diff = pright - pleft;
    MatrixXd prnew = pleft + ratio * diff;
    cout << "prnew" << prnew << endl;
    
    // tl.p3.resize(3,2);
    tl.p3 = Hcat(pleft, prnew);
    tl.r1 = 0.065;
    tl.r2 = 0.05;
    tl.r3 = 0.03;
    } 
};


#endif
