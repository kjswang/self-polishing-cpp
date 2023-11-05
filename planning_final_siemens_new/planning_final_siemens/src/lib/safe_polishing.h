#ifndef _SAFEPOLISHING_H_
#define _SAFEPOLISHING_H_

/*
* Safe polishing function definition
* polishing libraries 
* Author -
* Date - 
*/

#include "customized_print.h"
#include "load_parameter.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include "robot_property_eigen.h"
#include "distance_constraint_3d.h"
#include "structure.h"
#include "numerical_gradient.h"
#include "CapPos.h"
#include <vector>
#include "math_pack.h"
#include <cmath>
#include "QuadProg++.hh"
#include "global_var.h"
// #include "ProblemCFS_qp3d_position.h"
// #include "PlanningProblem_3d_position.h"
using namespace std;




void addLowUpCons(MatrixXd& A, MatrixXd& b, MatrixXd& x, MatrixXd lb, MatrixXd ub){
    /* 
    * this is for quadratic programming 
    * add the lower bound and upper bound for variable x 
    * lb: is the lower bound 
    * up: is the upper bound
    */
    int dim = x.rows()+2;
    // set the lower bound
    MatrixXd Alow = -1*MatrixXd::Identity(dim,dim);

    assert(lb.cols() == 1); // make sure bl is a column vector
    MatrixXd blow = -1*lb;
    // set the upper bound
    MatrixXd Aup = MatrixXd::Identity(dim,dim);

    assert(ub.cols() == 1);
    MatrixXd bup = ub;
    // concatenate the upper bound condition and lower bound condition together with A, b
    // concatenate with lower bound
    A = Vcat(A,Alow);
    b = Vcat(b,blow);
    // concatenate with upper bound
    A = Vcat(A,Aup);
    b = Vcat(b,bup);
    // all is done now .
}


MatrixXd Jacobi(MatrixXd theta, MatrixXd DH, int n, MatrixXd p){
    /* the function calculates the Jacobi matrix (J) of
    * a given point p (in the base coordinate) on the link n
    */
    MatrixXd z, r_0, J, TCP_T;
    int nlink = theta.rows();
    DH.block(0,0,nlink,1) = theta;
    z = MatrixXd::Zero(3,n); // z(:,i) is the z axis of link i in the base coordinate
    r_0 = MatrixXd::Zero(3,n+1); // r_0(:,i) is the coordinate of the i-th origin
    J = MatrixXd::Zero(nlink,n); // Jacobian
    TCP_T = MatrixXd::Identity(4,4);

    MatrixXd alpha, A, D, q;
    alpha = DH.col(3);
    A = DH.col(2);
    D = DH.col(1);
    q = DH.col(0);

    MatrixXd tmp;
    for (int i=0; i<n; ++i){
        z.col(i) = TCP_T.block(0,2,3,1);
        r_0.col(i) = TCP_T.block(0,3,3,1);
        tmp.resize(4,4);
        tmp << cos(q(i,0)), -sin(q(i,0))*cos(alpha(i,0)), sin(q(i,0))*sin(alpha(i,0)), A(i,0)*cos(q(i,0)),
                sin(q(i,0)), cos(q(i,0))*cos(alpha(i,0)), -cos(q(i,0))*sin(alpha(i,0)), A(i,0)*sin(q(i,0)),
                0, sin(alpha(i,0)), cos(alpha(i,0)), D(i,0),
                0, 0, 0, 1;
        TCP_T = TCP_T*tmp;
    }
    r_0.col(n) = TCP_T.block(0,3,3,1);
    // TCP_T=TCP_T*... % up to now, forward kinematics, S0 to S6
    MatrixXd up;
    for (int i=0; i<n; ++i){
        up = Matcross(r_0.col(i) - p, z.col(i));
        J.col(i) = Vcat(up,z.col(i));
    }

    return J;
}


MatrixXd readXYZFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return MatrixXd();
    }

    std::string line;
    std::vector<Vector3d> points;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
            std::cerr << "Invalid data format in file: " << filename << std::endl;
            return MatrixXd();
        }
        points.push_back(Vector3d(x, y, z));
    }

    file.close();

    if (points.empty()) {
        std::cerr << "No points found in file: " << filename << std::endl;
        return MatrixXd();
    }

    MatrixXd pointMatrix(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        pointMatrix.row(i) = points[i].transpose();
    }

    return pointMatrix;
}


void load_PC(MatrixXd& PC, MatrixXd& PC_idx){
    string sample_file = "parameter/sampled_wp.xyz";
    string xyz_file = "parameter/PC_idx.xyz";
    PC = readXYZFile(sample_file);
    load_PC_idx(PC_idx);
}





MatrixXd ForKine(MatrixXd theta_ini, MatrixXd DH, MatrixXd base, capsule RoCap[]){
    /*
    * Brief version of forward kinematics 
    * only give the end effector position of robot
    * the last joint coordinate frame position is considered as end effector position
    */
    MatrixXd Msix2tool;
    Msix2tool.resize(4,4);
    if (planner == 1){ // polishing planner
        loadjnt2tool(Msix2tool); // this is measured transformation matrix from robot end effector frame to tool tip frame
    }
    if (planner == 0){ // measurement planner
        loadjnt2laser(Msix2tool); // this is measured transformation matrix from robot end effector frame to laser frame
    }
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
    return epos;
}

double dist_nullspace(MatrixXd x0, MatrixXd theta0, MatrixXd nullspace, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, double alpha){
    MatrixXd theta;
    MatrixXd zero_tmp(1,1);
    zero_tmp(0,0) = 0;

    theta = theta0 + alpha*(nullspace*x0);

    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;
    lineseg* pos = CapPos_real(base, DH, cap, tl);

    MatrixXd point1, point2;
    double dist,radius;
    vector<double> distList(sizeof(pos));
    
    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC(point1, point2, PC.transpose());
        dist -= radius;

        distList[i] = dist;
    }

    MatrixXd c_end = ForKine(theta, DH, base, cap);
    MatrixXd c_last = pos[7].p2;
    MatrixXd v = c_end - c_last;
    double norm = v.norm();
    MatrixXd normV3d = v/norm;
    MatrixXd c_end1 = c_end - 0.03 * normV3d;

    double dist_end = distLinSeg2PC(c_end1, c_end1, PC.transpose()) - 0.03;
    distList[8] = dist_end;

    auto startIterator = distList.begin() + 3; // Points to the 4th element
    auto endIterator = distList.end();         // Points to the end

    d = *std::min_element(startIterator, endIterator);
    return -d;
}

void setRobotGoal(Robot& robot, int& nstep1, int& nstep2, int& nwait, int& tfinal, MatrixXd& theta_ini, MatrixXd& weldTraj, MatrixXd& center_point){
    /*
    * Set the reference end effector trajectories for robot
    * start2exe: the steps from initial pos to executation pos 
    * nwait: the steps to wait before execution 
    * tfinal: the reference weld point number 
    */
    // make sure weldTraj dimension is correct
    assert(weldTraj.rows() == 3);

    // get pos initial 

    MatrixXd posini = ForKine(theta_ini, robot.DH, robot.base, robot.cap);
    MatrixXd posexe = weldTraj.block(0,0,3,1);

    // set robot end effector refrence trajectory
    // from the initial states to execution states 
    MatrixXd diff;
    MatrixXd diff1 = center_point - posini;
    MatrixXd diff2 = posexe - center_point;
    MatrixXd tmp;
    for (int t=1; t<=nstep1; ++t){
        tmp = posini + double(t)/nstep1*diff1;
        robot.goal = Hcat(robot.goal, tmp);
    }

    for (int t=1; t<=nstep2; ++t){
        tmp = center_point + double(t)/nstep2*diff2;
        robot.goal = Hcat(robot.goal, tmp);
    }

    // wait until the velocity and acceleration approaches 0
    for (int t=1; t<=nwait; ++t){
        robot.goal = Hcat(robot.goal, posexe);
    }

    // execute polishing (end-effector located precisely on the weld points)
    int ministep = 1; // interpolation (augment weld points)
    for (int t=0; t<tfinal-1; ++t){
        diff = weldTraj.block(0,t+1,3,1) - weldTraj.block(0,t,3,1);
        // interplation
        for (int st=1; st<=ministep; ++st){
            tmp = weldTraj.block(0,t,3,1) + double(st)/ministep*diff;
            robot.goal = Hcat(robot.goal, tmp);
        }
    }
    tfinal = (tfinal - 1)*ministep; // update tfinal steps

}

tuple<MatrixXd, MatrixXd> safetrack_auto_robot(MatrixXd theta0, MatrixXd wp_pos0, Robot robot, MatrixXd c1, MatrixXd PC, MatrixXd PC_idx){
    MatrixXd PC_new, wp_pos_new, c1_new, c0;
    tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos0);
    wp_pos_new = wp_pos0;
    c1_new = next_point_WP(wp_pos0, c1, PC);
    c0 = ForKine(theta0, robot.DH, robot.base, robot.cap);

    assert(robot.nlink == 6);
    double penalty1[6] = {5,5,2,1,5,10};
    double penalty2[6] = {5,5,100,100,100,100};
    double dist;

    MatrixXd H = MatrixXd::Zero(robot.nlink, robot.nlink);
    if (PC_idx.isZero()){
        for (int jnt=0; jnt<robot.nlink; ++jnt){
        H(jnt,jnt) = penalty2[jnt];
        }
    }
    else {
        for (int jnt=0; jnt<robot.nlink; ++jnt){
        H(jnt,jnt) = penalty1[jnt];
        }
    }

    MatrixXd f = -H.transpose()*theta0;

    MatrixXd Sstack, Lstack, ref_grad, l, s, distMat;
    MatrixXd abs_PC_idx = PC_idx.array().abs();
    int max_value = abs_PC_idx.maxCoeff();

    for (int j = 1; j <= max_value; j++) {
        MatrixXd curPC;
        for (int i = 0; i < PC_idx.rows(); i++) {
            if (PC_idx(i) == j){
                curPC.conservativeResize(curPC.rows() + 1, curPC.cols());
                curPC.row(curPC.rows() - 1) = PC_new.row(i);
            }
        }

        dist = dist_arm_PC(theta0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
        
        distMat.resize(1,1);
        distMat << dist;
        ref_grad = central_diff_polish(theta0, robot.DH, robot.base, robot.cap, robot.tl, curPC);

        s = distMat - ref_grad*theta0;
        l = -ref_grad;
        Sstack = Vcat(Sstack, s);
        Lstack = Vcat(Lstack, l);
    }

    MatrixXd Jac, Diff;
    Jac = Jacobi(theta0, robot.DH, robot.nlink, c0);
    Diff = Jac.block(0,0,3,Jac.cols());

    MatrixXd Aeq, beq;
    Aeq = Diff;
    beq = c1_new - c0 + Diff*theta0;


        // solve QP 
    // set the QP objective and linear and quadratic term 
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;

    // set objective quadratic and linear terms
    setObjValue(G, g0, H, f); // G is quadratic term and f is the linear term 

    // set linear inequality constraints    
    if (Lstack.rows()==0){// when inequality is not applicable 
        CI.resize(robot.nlink,0);
        ci0.resize(0);
    }
    else{
        MatrixXd nLT = -1*Lstack.transpose(); // negative transpose
        setConstraint(CI, ci0, nLT, Sstack); // set constriants
    }

    // set linear equality cosntraints
    if (beq.rows()==0){// when equality is not applicable 
        CE.resize(robot.nlink,0);
        ce0.resize(0);
    }
    else{
        MatrixXd nAeqT = -1*Aeq.transpose(); // negative transpose
        setConstraint(CE, ce0, nAeqT, beq); // set cosntriants
    }
    // solve the QP problem 
    // set the x as initial point 
    // // check the correctness of variables
    // cout << "H is :" << H << endl;
    // cout << "f is :" << f << endl;
    // if (beq.rows()==0){// when equality is not applicable 
    //     cout << "CE is :" << 0 << endl;
    //     cout << "ce0 is :" << 0 << endl;
    // }
    // else{
    //     cout << "CE is :" << Aeq << endl;
    //     cout << "ce0 is :" << beq << endl;
    // }
    
    // if (Lstack.rows()==0){// when inequality is not applicable 
    //     cout << "CI is :" << 0 << endl;
    //     cout << "ci0 is :" << 0 << endl;
    // }
    // else{
    //     cout << "CI is :" << Lstack << endl;
    //     cout << "ci0 is :" << Sstack << endl;
    // }
    

    QPxset(x, theta0);// set initial value of u;
    double tmp_cost = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    cout << "the temporal cost is " << tmp_cost << endl;
    cout << "valid quadprog is solved" << endl;
    MatrixXd theta_new;
    quad2matrix(x, theta_new); // push the solved x value to theta_new
    // MatrixXd iediff = Sstack - Lstack*theta_new;
    // MatrixXd ediff = beq - Aeq*theta_new;
    // make sure the inequality and equality constriants are satisfied perfectly
    return std::make_tuple(theta_new, wp_pos_new);
}

bool check_ineq_need(MatrixXd theta, MatrixXd DH, MatrixXd base, capsule cap[], double consider_line){
    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;
    lineseg pos[DH.rows()];
    MatrixXd M[DH.rows()+1];
    CapPos(base, DH, cap, M, pos);
    double line_x = consider_line;
    int cross_line_num = 0;


    for (int i = 0; i < nlink; i++) {
        if(pos[i].p1(0,0) < line_x && pos[i].p2(0,0) < line_x){
            continue;
        }
        else {
            cross_line_num++;
        }
    }
   
    
    return cross_line_num > 0;

}

tuple<MatrixXd, MatrixXd> safetrack_auto_collaboration(MatrixXd theta0, MatrixXd wp_pos0, Robot robot, MatrixXd c1, MatrixXd PC, MatrixXd PC_idx){
    MatrixXd PC_new, wp_pos_new, c1_new, c0, x0;
    //tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos0);
    //wp_pos_new = wp_pos0;
    c1_new = next_point_WP(wp_pos0, c1, PC);
    double consider_line = 1.2;
    c0 = ForKine(theta0, robot.DH, robot.base, robot.cap);

    x0 = Vcat(theta0,wp_pos0.row(1));
    x0 = Vcat(x0,wp_pos0.row(2));

    assert(robot.nlink == 6);
    double penalty[8] = {5,5,2,1,5,10,50,40};
    double dist;

    MatrixXd H = MatrixXd::Zero(8, 8);
    for (int jnt=0; jnt<8; ++jnt){
        H(jnt,jnt) = penalty[jnt];
        }
 
    MatrixXd f = -H.transpose()*x0;


    MatrixXd Sstack, Lstack, ref_grad, l, s, distMat;
    MatrixXd abs_PC_idx = PC_idx.array().abs();
    int max_value = abs_PC_idx.maxCoeff();
    bool need_flag = check_ineq_need(theta0, robot.DH, robot.base, robot.cap, consider_line);
    if (need_flag){
    for (int j = 1; j <= max_value; j++) {
        MatrixXd curPC;
        for (int i = 0; i < PC_idx.rows(); i++) {
            if (PC_idx(i) == j){
                curPC.conservativeResize(curPC.rows() + 1, curPC.cols());
                curPC.row(curPC.rows() - 1) = PC.row(i);
            }
        }

        dist = dist_arm_PC_WP(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
        
        distMat.resize(1,1);
        distMat << dist;
        ref_grad = central_diff_polish_collaboration(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);

        s = distMat - ref_grad*theta0;
        l = -ref_grad;
        Sstack = Vcat(Sstack, s);
        Lstack = Vcat(Lstack, l);
    }
    }
    MatrixXd Jac, Diff_Jac, Diff_cfun, Diff;
    Jac = Jacobi(theta0, robot.DH, robot.nlink, c0);

    Diff_Jac = Jac.block(0,0,3,Jac.cols());
    Diff_cfun = Diff_Jac_num_grad(c1, PC, wp_pos0);
    Diff_cfun = Diff_cfun.block(0, 1, Diff_cfun.rows(), 2);

    Diff = Hcat(Diff_Jac,Diff_cfun);

    MatrixXd Aeq, beq;
    Aeq = Diff;
    beq = c1_new - c0 + Diff*x0;


    MatrixXd lb = MatrixXd::Zero(8,1);
    MatrixXd ub = MatrixXd::Zero(8,1);
    lb.block(0,0,robot.lb.rows(),1) << robot.lb;
    lb(6,0) = -1.0;
    lb(7,0) = -1.0;

    ub.block(0,0,robot.ub.rows(),1) << robot.ub;
    ub(6,0) = 1.0;
    ub(7,0) = 1.0;

    // MatrixXd lb = MatrixXd::Zero(6,1);
    // MatrixXd ub = MatrixXd::Zero(6,1);
    // lb.block(0,0,robot.lb.rows(),1) << robot.lb;

    // ub.block(0,0,robot.ub.rows(),1) << robot.ub;

    addLowUpCons(Lstack, Sstack, theta0, lb, ub);

    //addLowUpCons(Lstack, Sstack, theta0, robot.lb, robot.ub);

        // solve QP 
    // set the QP objective and linear and quadratic term 
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;

    // set objective quadratic and linear terms
    setObjValue(G, g0, H, f); // G is quadratic term and f is the linear term 

    // set linear inequality constraints    
    if (Lstack.rows()==0){// when inequality is not applicable 
        CI.resize(robot.nlink,0);
        ci0.resize(0);
    }
    else{
        MatrixXd nLT = -1*Lstack.transpose(); // negative transpose
        setConstraint(CI, ci0, nLT, Sstack); // set constriants
    }

    // set linear equality cosntraints
    if (beq.rows()==0){// when equality is not applicable 
        CE.resize(robot.nlink,0);
        ce0.resize(0);
    }
    else{
        MatrixXd nAeqT = -1*Aeq.transpose(); // negative transpose
        setConstraint(CE, ce0, nAeqT, beq); // set cosntriants
    }
    // solve the QP problem 
    // set the x as initial point 
    // // check the correctness of variables
    // cout << "H is :" << H << endl;
    // cout << "f is :" << f << endl;
    // if (beq.rows()==0){// when equality is not applicable 
    //     cout << "CE is :" << 0 << endl;
    //     cout << "ce0 is :" << 0 << endl;
    // }
    // else{
    //     cout << "CE is :" << Aeq << endl;
    //     cout << "ce0 is :" << beq << endl;
    // }
    
    // if (Lstack.rows()==0){// when inequality is not applicable 
    //     cout << "CI is :" << 0 << endl;
    //     cout << "ci0 is :" << 0 << endl;
    // }
    // else{
    //     cout << "CI is :" << Lstack << endl;
    //     cout << "ci0 is :" << Sstack << endl;
    // }
    
    // cout << G << endl;
    // cout << g0 << endl;
    // cout << CE << endl;
    // cout << ce0 << endl;
    // cout << CI << endl;
    // cout << ci0 << endl;
    // cout << x << endl;
    QPxset(x, theta0);// set initial value of u;
    // cout << G << endl;
    // cout << g0 << endl;
    // cout << CE << endl;
    // cout << ce0 << endl;
    // cout << CI << endl;
    // cout << ci0 << endl;
    // cout << x << endl;

    double tmp_cost = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    cout << "the temporal cost is " << tmp_cost << endl;
    cout << "valid quadprog is solved" << endl;
    MatrixXd x_new, theta_new;
    quad2matrix(x, x_new); // push the solved x value to theta_new
    theta_new = x_new.topRows(6);
    wp_pos_new = x_new.bottomRows(2);

    MatrixXd wp_pos_new2(3,1);
    wp_pos_new2 << 0,wp_pos_new;
    // MatrixXd iediff = Sstack - Lstack*theta_new;
    // MatrixXd ediff = beq - Aeq*theta_new;
    // make sure the inequality and equality constriants are satisfied perfectly
    return std::make_tuple(theta_new, wp_pos_new2);
}

bool check_collision_complete_PC_cluster(MatrixXd theta0, Robot robot, MatrixXd PC, MatrixXd PC_idx){
    bool col_flag = false;
    MatrixXd abs_PC_idx = PC_idx.array().abs();
    int max_value = abs_PC_idx.maxCoeff();
    double minDist = INFINITY;
    double dist;
    for (int j=1; j<=max_value;j++){

        vector<int> indices;
        for (int i = 0; i < PC_idx.cols(); ++i){
            if(PC_idx(0,i) == j){
                indices.push_back(i);
            }
        }

        MatrixXd curPC(indices.size(),3);
        for (int i =0; i < indices.size(); i++){
            curPC.row(i) = PC.row(indices[i]);
        }

        dist = dist_arm_PC(theta0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
        minDist = min(dist, minDist);
        if (minDist <= 0){
            col_flag = true;
            cout << "collision happens" << endl;
            break;
        }
    }
    return col_flag;
}

bool checkFeasible(MatrixXd theta, MatrixXd PC, MatrixXd PC_idx, MatrixXd DH, MatrixXd base, capsule RoCap[], tool tl){
    double PC_minX = PC.col(0).minCoeff();
    double PC_maxX = PC.col(0).maxCoeff();
    double PC_minY = PC.col(1).minCoeff();
    double PC_maxY = PC.col(1).maxCoeff();
    double PC_minZ = PC.col(2).minCoeff();
    double PC_maxZ = PC.col(2).maxCoeff();
    bool feasible = true;

    double offset = 0.2;
    // MatrixXd PC_hull = PC.rowwise().any() && (PC_idx.array() < 0).colwise().all();
    // MatrixXd PC_hull = (PC.rowwise().any()).cast<double>() && (PC_idx.array() < 0).colwise().all().cast<double>();

    vector<int> indices;
    for (int i = 0; i < PC_idx.cols(); ++i){
        if(PC_idx(0,i) < 0){
            indices.push_back(i);
        }
    }

    MatrixXd PC_hull(indices.size(),3);
    for (int i =0; i < indices.size(); i++){
        PC_hull.row(i) = PC.row(indices[i]);
    }

    // MatrixXd PC_hull;
    // VectorXd idx_mask = (PC_idx.array() < 0).cast<double>();
    // PC_hull = PC.array().rowwise() * idx_mask.transpose().array();

    DH.col(0) = theta;
    //MatrixXd pos, M;
    // lineseg pos[DH.rows()];
    // MatrixXd M[DH.rows()+1];
    // CapPos(base,DH,RoCap,M,pos);
    lineseg* pos = CapPos_real(base, DH, RoCap, tl);

    for (int i = 0; i < sizeof(pos); i++) {
        MatrixXd point1 = pos[i].p1;
        MatrixXd point2 = pos[i].p2;
        double radius = pos[i].r;
        if (point1(0) > (PC_minX + offset) && inhull(point1.transpose(), PC_hull) == 0) {
            feasible = false;
            return feasible;
        }
        if (point2(0) > (PC_minX + offset) && inhull(point2.transpose(), PC_hull) == 0) {
            feasible = false;
            return feasible;
        }
    }

    MatrixXd c_end = ForKine(theta, DH, base, RoCap);
    MatrixXd c_last = pos[7].p2;
    MatrixXd v = c_end - c_last;
    double norm = v.norm(); 
    MatrixXd normV3d = v / norm;
    MatrixXd c_end1 = c_end - 0.03 * normV3d;
    if (distLinSeg2PC(c_end1, c_last, PC.transpose()) < 0.03 || distLinSeg2PC(c_end, c_last, PC.transpose()) < 0) {
        feasible = false;
        return feasible;
    }

    return feasible;
}

bool checkValid(MatrixXd theta_0, MatrixXd theta_1, MatrixXd wp_pos_0, MatrixXd wp_pos_1, MatrixXd c1, Robot robot, MatrixXd PC_origin, MatrixXd PC_idx){
    bool valid = true;
    MatrixXd c0 = ForKine(theta_1, robot.DH, robot.base, robot.cap);
    MatrixXd PC_0, PC_1;
    tie(PC_0, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_0);
    tie(PC_1, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_1);
    MatrixXd diff = c1 - c0;
    if (diff.norm() > 0.001){
        cout << "P0" << endl;
        valid = false;
    }
    cout << check_collision_complete_PC_cluster(theta_1,robot,PC_1,PC_idx) << endl;
    if (check_collision_complete_PC_cluster(theta_1,robot,PC_1,PC_idx) == true){
        cout << "P1" << endl;
        valid = false;
    }

    if (checkFeasible(theta_1, PC_1, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0){
        cout << "P2" << endl;
        valid = false;
    }

    MatrixXd diff_wp_pos = wp_pos_1 - wp_pos_0;
    MatrixXd diff_theta = theta_1 - theta_0;
    int steps = 20;
    double minDist = INFINITY;

    MatrixXd PC_step = PC_0;
    for (int s = 0; s <= steps; s++){
        MatrixXd theta_step = theta_0 + s * diff_theta / steps;
        if (diff_wp_pos.norm() > 0) {
            MatrixXd wp_pos_step = wp_pos_0 + s * diff_wp_pos/steps;
            tie(PC_step, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_step);
        }

        bool col_flag = check_collision_complete_PC_cluster(theta_step, robot, PC_step, PC_idx);
        if (col_flag == 1) {
            valid = 0;
            cout << "P3" << endl;
            break;
        }

        if (checkFeasible(theta_step, PC_step, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0) {
            valid = 0;
            cout << "P4" << endl;
            break;
        }

    }
    return valid;
}

MatrixXd iterTrack(MatrixXd theta0, MatrixXd c1, Robot robot){
    for (int i = 0; i < 100; i++){
        MatrixXd c0 = ForKine(theta0, robot.DH, robot.base, robot.cap);
        MatrixXd dc = c1-c0;
        if(dc.norm() < 0.001){
            break;
        }

        MatrixXd Jac = Jacobi(theta0, robot.DH, robot.nlink, c0-robot.base);
        cout << Jac << endl;
        MatrixXd Jac_new = Jac.topRows(3);
        MatrixXd dtheta = Jac_new.completeOrthogonalDecomposition().pseudoInverse() * dc;
        theta0 = theta0 + dtheta;
    }
    return theta0;
}

template <typename Number> // 'Number' can be 'double' or 'std::complex<double>'
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> kernel_COD(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  Eigen::CompleteOrthogonalDecomposition<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      cod;
  cod.compute(M);
  unsigned rk = cod.rank();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> P =
      cod.colsPermutation();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> V =
      cod.matrixZ().transpose();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Kernel =
      P * V.block(0, rk, V.rows(), V.cols() - rk);
  return Kernel;
}



void proc_nullspace_opt(MatrixXd theta_cur, MatrixXd wp_pos_cur, Robot robot, MatrixXd PC_origin, MatrixXd PC_idx, MatrixXd c_next, int safe_planner, string nullspace_mode){
    double alpha0 = 0.1;
    MatrixXd theta_tmp = theta_cur;
    MatrixXd wp_pos_tmp = wp_pos_cur;
    MatrixXd x0(3,1);
    x0 << 0,0,0;
    MatrixXd c0 = ForKine(theta_tmp, robot.DH, robot.base, robot.cap);
    MatrixXd PC_0;
    tie(PC_0, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_cur);

    MatrixXd c1 = next_point_WP(wp_pos_cur, c_next, PC_origin);
    theta_tmp = iterTrack(theta_cur, c1, robot);
    MatrixXd Jac = Jacobi(theta_tmp, robot.DH, robot.nlink, c0 - robot.base);
    MatrixXd Jac_new = Jac.topRows(3);
    MatrixXd Jac_new_T = Jac_new.transpose();
    MatrixXd Jac_null = kernel_COD(Jac_new);

    // FullPivLU<MatrixXd> lu(Jac_new);
    // MatrixXd A_null_space = lu.kernel();
    // cout << A_null_space << endl;

    // HouseholderQR<MatrixXd> qr(A_null_space);
    // MatrixXd orthonormalBasis = qr.householderQ();
    // cout << orthonormalBasis << endl;


    // JacobiSVD<MatrixXd> svd(Jac_new, ComputeFullU | ComputeFullV);
    // MatrixXd U = svd.matrixU();
    // MatrixXd S = svd.singularValues().asDiagonal();
    // MatrixXd V = svd.matrixV();

    // cout << U << endl;
    // cout << S << endl;
    // cout << V << endl;
    // const double epsilon = 1e-9;
    // int nullCols = 0;
    // for (int i = 0; i < svd.singularValues().size(); ++i){
    //      if(svd.singularValues()(i) > epsilon){
    //          nullCols++;
    //      }
    // }

    // MatrixXd Jac_null = svd.matrixV().rightCols(6 - nullCols);
    // // cout << svd.matrixV() << endl;
    // // cout << Jac_null << endl;
    // Jac_null.transposeInPlace();
    // cout << Jac_null << endl;
    // HouseholderQR<MatrixXd> qr(Jac_null.transpose());
    // MatrixXd orthonormalBasis = qr.householderQ().transpose();
    // std::cout << "Orthonormal basis for null space:\n" << orthonormalBasis << std::endl;

    // FullPivLU<MatrixXd> lu(Jac_new);
    // MatrixXd Jac_null = lu.kernel();
    // cout << "Jac_null" << Jac_null << endl;
    // MatrixXd orthonormalBasis = Jac_null.householderQr().householderQ();
    // cout << orthonormalBasis << endl;


    // FullPivLU<MatrixXd> lu(Jac);
    // MatrixXd Jac_null = lu.kernel();

    double dist = dist_nullspace(x0, theta_tmp, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
    cout << dist << endl;

    vector<double> d_list;
    d_list.push_back(dist);
    double d_min = dist;
    double alpha = alpha0;
    MatrixXd theta_opt;

    for (int i=1; i<=10; i++) {
        if(dist < 0 && safe_planner ==0 && checkValid(theta_cur, theta_tmp, wp_pos_cur, wp_pos_cur, c1, robot, PC_origin, PC_idx) == 1){
            theta_opt = theta_tmp;
            break;
        }
        c0 = ForKine(theta_tmp, robot.DH, robot.base, robot.cap);
        Jac = Jacobi(theta_tmp, robot.DH, robot.nlink, c0-robot.base);
        Jac_new = Jac.topRows(3);
        Jac_null = kernel_COD(Jac_new);
        if (nullspace_mode == "cmaes"){
            continue;
        }
        else{
            MatrixXd A, b, Aeq, beq, lb, ub;


        }

    }

}

void safe_polish_procedure(MatrixXd& safe_theta, MatrixXd& safe_traj, MatrixXd& safe_wp_pos, Robot robot, int nstep1, int nstep2, int nwait, int tfinal, 
                MatrixXd &theta_pre, MatrixXd &c_pre, double thres, MatrixXd &M_PC_0, MatrixXd &wp_pos_pre, string track_mode_main, 
                string track_solver, MatrixXd PC_origin, MatrixXd PC_idx, string explore_mode, string nullspace_mode){

MatrixXd c_next, diff;
bool col_flag;
MatrixXd theta_new, wp_pos_new;
int repeat_cnt;

vector<double> curLog(5);
vector<vector<double>> log;

for (int t=0; t<tfinal + nstep1 + nstep2 + nwait; ++t){
    cout << "---------- time step " << t << " ----------" << endl;
    c_next = robot.goal.block(0,t,3,1);

    MatrixXd inv = M_PC_0.inverse();
    c_next = setVertice(c_next.transpose(), inv).transpose();

    curLog[0] = t;

    col_flag = true;

    MatrixXd theta_cur = theta_pre;
    MatrixXd c_cur = c_pre;
    MatrixXd wp_pos_cur = wp_pos_pre;
    MatrixXd PC_new;

    int iteration_cnt = 1;

    while (iteration_cnt <= 20)
    {
        iteration_cnt++;

        if (track_mode_main == "robot")
        {
            if (track_solver == "ICOP")
            {
                cout << "****************Robot ICOP Tracking*****************" << endl;
                 tie(theta_new, wp_pos_new) = safetrack_auto_robot(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
            }
            // else 
            // {
            //     cout << "****************Robot Fmincon Tracking*****************" << endl;
            //     tie(theta_new, wp_pos_new) = safetrack_auto_fmincon(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
            // }
        }

        if (track_mode_main == "collaboration")
        {
            cout << "****************Collabration Tracking*****************" << endl;
            tie(theta_new, wp_pos_new) = safetrack_auto_collaboration(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
        }       

        if (theta_new.size() == 0)
            continue;

        MatrixXd c_new = ForKine(theta_new, robot.DH, robot.base, robot.cap);
        tie(PC_new, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_new);
        MatrixXd c_target = next_point_WP(wp_pos_new, c_next, PC_origin);

        c_pre = c_new;
        theta_pre = theta_new;
        wp_pos_pre = wp_pos_new;

        col_flag = check_collision_complete_PC_cluster(theta_new, robot, PC_new, PC_idx);
        diff = c_target - c_pre;
        if (diff.norm() <= thres && col_flag == 0){
            break;
        }
    }

    curLog[1] = iteration_cnt;
    int need_sample = 0;
    int stop_flag = 0;
    int success = 0;
    int safe_planner = 0;

    if ( t < 10-1 && t > 5-1) {
        safe_planner = 1;
    }

    
    MatrixXd c1 = next_point_WP(wp_pos_pre, c_next, PC_origin);
    if (checkValid(theta_cur, theta_pre, wp_pos_cur, wp_pos_pre, c1, robot, PC_origin, PC_idx) == 0){
        safe_planner = 1;
    }

    if (safe_planner || iteration_cnt > 20 || checkFeasible(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0){
        if (explore_mode == "None") {
            stop_flag = 1;
        }
        else {
            need_sample = 1;
            success = 0;

            // if (explore_mode == "CSsample") {
            //     proc_auto_sample;
            // }
            // if (explore_mode == "NSsample") {
            //     proc_auto_sample_nullspace;
            // }
            if (explore_mode == "NSopt") {
                proc_nullspace_opt(theta_cur, wp_pos_cur, robot, PC_origin, PC_idx, c_next, safe_planner, nullspace_mode);
                continue;
            }
            
            if (success == 0){
                stop_flag =1;
            }
            else {
                tie(PC_new, ignore, ignore, ignore) = processPC(PC_origin, wp_pos_pre);
                if (checkFeasible(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap,robot.tl) == 0) {
                    stop_flag = 1;
                }
            }
        }
    }

    curLog[2] = 0; //T1
    curLog[3] = 0; //T2
    curLog[4] = need_sample;
    log.push_back(curLog);

    if (stop_flag == 1) {
            break;
        }
    safe_traj = Hcat(safe_traj, c_pre);
    safe_theta = Hcat(safe_theta, theta_pre);
    safe_wp_pos = Hcat(safe_wp_pos, wp_pos_pre);

    cout << "safe_traj" << safe_traj <<endl;
    cout << "safe_theta" << safe_theta << endl;
    cout << "safe_wp_pos" << safe_wp_pos << endl;
    }    
}





void safe_polish(MatrixXd &start2exe_traj, MatrixXd &execution_traj, MatrixXd &exit_traj){

string name = "GP50";
Robot robot(name);
cout << "DH:\n" << robot.DH << endl;
cout << "nliink:\n" << robot.nlink << endl;


// load_PC
MatrixXd PC, PC_idx;
load_PC(PC, PC_idx);



MatrixXd arr_axis3;
loadWeldTraj(arr_axis3);
MatrixXd arr_axis3T = arr_axis3.transpose();

MatrixXd abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top;
loadWorkpieceSetting(abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top);



MatrixXd weld_in_T = weld_in.transpose();
MatrixXd in_center = calculateRowMeans(weld_in_T);
MatrixXd radius = weld_in_T.colwise() - in_center.col(0);
VectorXd vecnorm = radius.colwise().norm();
radius = 0.035*radius;
radius.transposeInPlace();
MatrixXd temp = radius.array().colwise()/vecnorm.array();
weld_in = weld_in - temp;
MatrixXd temp1(1,3);
temp1 << 0,0.02,0;
MatrixXd temp2(1,3);
temp2 << 0,0,0.02;
weld_right = weld_right.rowwise() + temp1.row(0);
weld_left = weld_left.rowwise() - temp1.row(0);
weld_bottom = weld_bottom.rowwise() + temp2.row(0);
weld_top = weld_top.rowwise() - temp2.row(0);
arr_axis3 = weld_right.transpose();




string exp_mode = "robot";
string track_mode_main = "collaboration";
string track_solver = "ICOP";
string track_mode_sample = "robot";
string explore_mode = "NSopt";
string nullspace_mode = "sqp";

vector<string> exp_mode_list = {"robot", "workpiece"};
vector<string> track_solver_list = {"ICOP"};
vector<string> explore_mode_list = {"None", "CSsample", "NSample", "NSopt"};
vector<string> nullspace_mode_list = {"sqp", "interior-point", "active-set", "cmaes"};

vector<string> config_list;
// VectorXd theta_init_default(6);
// theta_init_default << 0, -2.3, 0.5, 0, 0.5, 0.0;
// Vector3d wp_pos_init_default(0.0, 0.4, 0.1);
// double alphaY_limit = 0.3; 
// double alphaZ_limit = 0.3;
// int nstep1 = 5;
// int nstep2 = 5;
// int nwait = 1;

double alphaY_limit, alphaZ_limit, thres;
int nstep1, nstep2, nwait;
MatrixXd theta_init_default(6,1), wp_pos_init_default(3,1);
loadSafePolishingSetting(alphaY_limit, alphaZ_limit, thres, nstep1, nstep2, nwait, theta_init_default, wp_pos_init_default);








vector<Vector3d> wp_pos_list;

double step = 0.1;

for (double alphaY = -alphaY_limit; alphaY <= alphaY_limit; alphaY += step) {
    for (double alphaZ = -alphaZ_limit; alphaZ <= alphaZ_limit; alphaZ += step){
            Vector3d wp_pos_cur(0.0, alphaY + 0.5, alphaZ);
            wp_pos_list.push_back(wp_pos_cur);
            MatrixXd PC_processed, M_PC_processed;
            MatrixXd base_point_processed, center_point_processed;
            tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_cur);
            MatrixXd arr_cur = setVertice(arr_axis3.transpose(), M_PC_processed).transpose();
    }

}

double joint2_limit = 0.2;
double joint3_limit = 0.2;
double joint4_limit = 0.2;
vector<VectorXd> theta_init_list;

for (double joint2 = -joint2_limit; joint2 <= joint2_limit; joint2 += step) {
    for (double joint3 = -joint3_limit; joint3 <= joint3_limit; joint3 += step) {
        for (double joint4 = -joint4_limit; joint4 <= joint4_limit; joint4 += step) {
            VectorXd theta_init(6);
            theta_init << 0.0, -2.3 + joint2, 0.5 + joint3, 0.0 + joint4, 0.5, -2.0;
            theta_init_list.push_back(theta_init);
            MatrixXd PC_processed, M_PC_processed;
            VectorXd base_point_processed, center_point_processed;
            tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_init_default);
            MatrixXd arr_cur = setVertice(arr_axis3.transpose(), M_PC_processed);
        }

    }

}

MatrixXd theta_init = theta_init_default;
MatrixXd wp_pos_init = wp_pos_init_default;
MatrixXd PC_processed, M_PC_processed;
MatrixXd base_point_processed, center_point_processed;
tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_init);
MatrixXd arr = setVertice(arr_axis3.transpose(), M_PC_processed).transpose();


int tfinal = arr.cols();

setRobotGoal(robot, nstep1, nstep2, nwait, tfinal, theta_init, arr, center_point_processed);

MatrixXd theta_pre, wp_pos_pre, c_pre, M_PC_0;
theta_pre = theta_init;
wp_pos_pre = wp_pos_init;
c_pre = ForKine(theta_init, robot.DH, robot.base, robot.cap);

M_PC_0 = M_PC_processed;

MatrixXd safe_traj = c_pre;
MatrixXd safe_theta = theta_pre;
MatrixXd safe_wp_pos = wp_pos_pre;

safe_polish_procedure(safe_theta, safe_traj, safe_wp_pos, robot, 
nstep1, nstep2, nwait, tfinal, theta_pre, c_pre, thres, M_PC_0, 
wp_pos_pre, track_mode_main, track_solver, PC, PC_idx, explore_mode, nullspace_mode);

MatrixXd safe_theta_T = safe_theta.transpose();
vector<VectorXd> uniqueRows;

for (int i = 0; i < safe_theta_T.rows(); ++i) {
    VectorXd currentRow = safe_theta_T.row(i);

    bool isUnique = true;
    for (const auto& row : uniqueRows) {
      if (currentRow.isApprox(row)) {
        isUnique = false;
        break;
      }
    }
    if (isUnique) {
      uniqueRows.push_back(currentRow);
    }
}

MatrixXd safe_theta_unique(uniqueRows.size(), safe_theta_T.cols());

for (int i = 0; i < uniqueRows.size(); ++i) {
    safe_theta_unique.row(i) = uniqueRows[i];
    }   

MatrixXd safe_theta_unique_final = safe_theta_unique.transpose();
}









#endif