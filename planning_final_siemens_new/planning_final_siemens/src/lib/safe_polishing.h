#ifndef _SAFEPOLISHING_H_
#define _SAFEPOLISHING_H_

/*
* Safe polishing function definition
* polishing libraries 
* Author -
* Date - 
*/

#include "customized_print.h"
#include "load_param.h"
#include <Eigen/Dense>
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

void safe_polish(MatrixXd &start2exe_traj, MatrixXd &execution_traj, MatrixXd &exit_traj){

string name = "GP50";
Robot robot(name);



// load_PC
MatrixXd PC;
load_PC(PC);

MatrixXd arr_axis3;
loadWeldTraj(arr_axis3);
MatrixXd arr_axis3T = arr_axis3.transpose();
MatrixXd abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top;
loadWorkpieceSettings(abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top)

MatrixXd weld_in_T = weld_in.transpose();
MatrixXd in_center = calculateMean(weld_in_T);
MatrixXd radius = weld_in_T.colwise() - in_center;
MatrixXd arr_axis3 = weld_right.transpose();

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
MatrixXd theta_init_default, wp_pos_init_default;
loadSafePolishingSetting(alphaY_limit, alphaZ_limit, nstep1, nstep2, theta_init_default, wp_pos_init_default);

vector<Vector3d> wp_pos_list;

double step = 0.1;

for (double alphaY = -alphaY_limit; alphaY <= alphaY_limit; alphaY += step) {
    for (double alphaZ = -alphaZ_limit; alphaZ <= alphaZ_limit; alphaZ += step){
            Vector3d wp_pos_cur(0.0, alphaY + 0.5, alphaZ);
            wp_pos_list.push_back(wp_pos_cur);
            MatrixXd PC_processed, M_PC_processed;
            VectorXd base_point_processed, center_point_processed;
            tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_cur);
            MatrixXd arr_cur = setVertice(arr_axis3T, M_PC);
            std::cout << "arr_cur:\n" << arr_cur << "\n";
    }

}

double joint2_limit = 0.2;
double joint3_limit = 0.2;
double joint4_limit = 0.2;
vector<Vector3d> theta_init_list = [];

for (double joint2 = -joint2_limit; joint2 <= joint2_limit; joint2 += step) {
    for (double joint3 = -joint3_limit; joint3 <= joint3_limit; joint3 += step) {
        for (double joint4 = -joint4_limit; joint4 <= joint4_limit; joint4 += step) {
            VectorXd theta_init(6);
            theta_init << 0.0, -2.3 + joint2, 0.5 + joint3, 0.0 + joint4, 0.5, -2.0;
            theta_init_list.push_back(theta_init);
            MatrixXd PC_processed, M_PC_processed;
            VectorXd base_point_processed, center_point_processed;
            tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_cur);
            MatrixXd arr_cur = setVertice(arr_axis3T, M_PC);
            std::cout << "arr_cur:\n" << arr_cur << "\n";
        }

    }

}

VectorXd theta_init = theta_init_default;
Vecotr3d wp_pos_init = wp_pos_init_default;
MatrixXd PC_processed, M_PC_processed;
VectorXd base_point_processed, center_point_processed;
tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_init);
MatrixXd arr = setVertice(arr_axis3.transpose(), M_PC)

int tfinal = arr.cols();
setRobotGoal(robot, nstep1, nstep2, nwait, tfinal, theta_init, arr, center_point_processed);

MatrixXd theta_pre, wp_pos_pre, c_pre, M_PC_0;
theta_pre = theta_init;
wp_pos_pre = wp_pos_init;
c_pre = ForKine(theta_init, robot.DH, robot.base, robot.cap);
M_PC_0 = M_PC;

MatrixXd safe_traj = c_pre;
MatrixXd safe_theta = theta_pre;
MatrixXd safe_wp_pos = wp_pos_pre;

safe_polish_procedure();

MatrixXd safe_theta_T = safe_theta.transpose();
vector<VectorXd> uniqueRows;

for (int i = 0; i < safe_theta_T.rows(); ++i) {
    VectorXd currentRow = transposed.row(i);

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


void load_PC(MatrixXd& PC){
    string sample_file = "parameter/sampled_wp.xyz";
    PC = readXYZFile(sample_file);
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

tuple<MatrixXd, MatrixXd, VectorXd, VectorXd> processPC(MatrixXd& PC, MatrixXd& PC_pos){
    MatrixXd M1;
    M1.resize(4,4);
    loadM1(M1);

    VectorXd x = PC.col(0);
    vector<int> idx;
    for (int i = 0; i < x.size(); ++i) {
        if (x(i) == 0) {
            idx.push_back(i);
        }
    }

    MatrixXd PC_x0(idx.size(), PC.cols());
    for (int i = 0; i < idx.size(); ++i) {
        PC_x0.row(i) = PC.row(idx[i]);
    }

    VectorXd center_point = calculateMean(PC_x0, 0);
    VectorXd mean_point = calculateMean(PC, 0);
    VectorXd base_point = mean_point;
    base_point(2) = PC.col(2).minCoeff() - 0.3;

    base_point = setVertice(base_point, M1);

    MatrixXd M_PC = M1;
    if (PC_pos.size() >= 3) {
        Vector3d elu = M1.block<3, 3>(0, 0).eulerAngles(0, 1, 2);
        for (int i = 0; i < 3; ++i) {
            Vector3d v = Vector3d::Zero();
            v(i) = 1;
            if (elu(i) == 0) {
                continue;
            }
            M_PC = applyRotation(M_PC, -v * elu(i), 1.0);
        }

        for (int i = 0; i < 3; ++i) {
            Vector3d v = Vector3d::Zero();
            v(i) = 1;
            if (PC_pos(i) == 0) {
                continue;
            }
            M_PC = applyRotation(M_PC, v * PC_pos(i), 1.0);
        }
    }

    PC = setVertice(PC, M_PC);
    center_point = setVertice(certer_point, M_PC);

    return std::make_tuple(PC.transpose(), M_PC, base_point, center_point);

}


MatrixXd setVertice(MatrixXd v, MatrixXd TransM) {
    MatrixXd result = v * TransM.block(0, 0, 3, 3).transpose();
    result.rowwise() += TransM.block(0, 3, 3, 1).transpose().replicate(v.rows, 1);
    return result;
}


void setRobotGoal(Robot& robot, int start2exe, int nwait, int& tfinal, MatrixXd& theta_ini, MatrixXd weldTraj, MatrixXd& center_point){
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


void safe_polish_procedure(MatrixXd& safe_theta, MatrixXd& safe_traj, Robot robot, MatrixXd planes, lineseg LineSegs[],  int start2exe, int nwait, int tfinal, 
                MatrixXd &theta_pre, MatrixXd theta0, MatrixXd &c_pre, double thres, MatrixXd anchor_point){

MatrixXd c_next, diff;
bool col_flag;
MatrixXd theta_new, wp_pos_new;
int repeat_cnt;

vector<double> curLog(5);
vector<vector<double>> log;

for (int t=0; t<tfinal + nwait + start2exe-1; ++t){
    cout << "---------- time step " << t << " ----------" << endl;
    c_next = robot.goal.block(0,t,3,1);
    c_next = setVertice(c_next.transpose(), M_PC_0.inverse());
    c_next = c_next.transpose();
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
            else 
            {
                cout << "****************Robot Fmincon Tracking*****************" << endl;
                tie(theta_new, wp_pos_new) = safetrack_auto_fmincon(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
            }
        }

        if (track_mode_main == "collaboration")
        {
            cout << "****************Collabration Tracking*****************" << endl;
            tie(theta_new, wp_pos_new) = safetrack_auto_collaboration(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
        }

        if (theta_new.size() == 0)
            continue;

        MatrixXd c_new = ForKine(theta_new, robotr.DH, robot.base, robot.cap);
        tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos_new);
        MatrixXd c_target = next_point_WP(wp_pos_new, c_next, PC);
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

    if ( t < 10 && t > 5) {
        safe_planner = 1;
    }

    MatrixXd c1 = next_point_WP(wp_pos_pre, c_next, PC);
    if (checkValid(theta_cur, theta_pre, wp_pos_cur, wp_pos_pre, c1, robot, PC_origin, PC_idx) == 0){
        safe_planner = 1;
    }

    if (safe_planner || iteration_cnt > 20 || checkFeasible(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap) == 0){
        if (explore_mode == "None") {
            stop_flag = 1;
        }
        else {
            need_sample = 1;
            success = 0;

            if (explore_mode == "CSsample") {
                proc_auto_sample;
            }
            if (explore_mode == "NSsample") {
                proc_auto_sample_nullspace;
            }
            if (explore_mode == "NSopt") {
                proc_nullspace_opt;
            }
            
            if (success == 0){
                stop_flag =1;
            }
            else {
                tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos_pre);
                if (checkFeasible(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap) == 0) {
                    stop_flag = 1;
                }
            }
        }
    }

    curLog(2) = T1;
    curLog(3) = T2;
    curLog(4) = need_sample;
    log.push_back(curLog);

    if (stop_flag == 1) {
            break;
        }
    safe_traj = Hcat(safe_traj, c_pre);
    safe_theta = Hcat(safe_theta, theta_pre);
    safe_wp_pos = Hcat(safe_wp_pos, wp_pos_pre);

    t++;
    }
    
    return 0;
}

MatrixXd next_point_WP(MatrixXd& wp_pos, MatrixXd& c1, MatrixXd& PC_origin) {
    MatrixXd curPC;
    MatrixXd M_PC;
    tie(curPC, M_PC, ignore, ignore) = processPC(PC_origin, wp_pos);
    MatrixXd c_next = setVertice(c1.transpose(), M_PC).transpose();
    return c_next;
}


tuple<MatrixXd, MatrixXd> safetrack_auto_robot(){
    MatrixXd PC_new, wp_pos_new, c1_new, c0;
    tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos0);
    wp_pos_new = wp_pos0;
    c1_new = next_point_WP(wp_pos0, c1, PC);
    c0 = ForKine(theta0, robot.DH, robot.base, robot.cap);

    assert(robot.nlink == 6);
    double penalty1[6] = {5,5,2,1,5,10};
    double penalty2[6] - {5,5,100,100,100,100};
    double dist;

    MatrixXd H = MatrixXd::Zero(robot.nlink, robot.nlink);
    if (PC_idx == 0){
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
        MaxtrixXd curPC;
        for (int i = 0; i < PC_idx.rows(); i++) {
            if (PC_idx(i) == j){
                curPC.convervativeResize(curPC.rows() + 1, curPC.cols());
                curPC.row(curPC.rows() - 1) = PC_new.row(i);
            }
        }

        dist = dist_arm_PC(theta0, robot.DH, robot.base, robot.cap, curPC);
        
        distMat.resize(1,1);
        distMat << dist;
        ref_grad = central_diff_polish(theta0, robot.DH, robot.base, robot.cap, curPC);

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


tuple<MatrixXd, MatrixXd> safetrack_auto_collaboration(){
    MatrixXd PC_new, wp_pos_new, c1_new, c0, x0;
    //tie(PC_new, ignore, ignore, ignore) = processPC(PC, wp_pos0);
    //wp_pos_new = wp_pos0;
    c1_new = next_point_WP(wp_pos0, c1, PC);
    c0 = ForKine(theta0, robot.DH, robot.base, robot.cap);

    x0 = Hcat(theta0,wp_pos0.row(1));
    x0 = Hcat(x0,wp_pos0.row(2));
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

    for (int j = 1; j <= max_value; j++) {
        MaxtrixXd curPC;
        for (int i = 0; i < PC_idx.rows(); i++) {
            if (PC_idx(i) == j){
                curPC.convervativeResize(curPC.rows() + 1, curPC.cols());
                curPC.row(curPC.rows() - 1) = PC.row(i);
            }
        }

        dist = dist_arm_PC_WP(x0, robot.DH, robot.base, robot.cap, curPC);
        
        distMat.resize(1,1);
        distMat << dist;
        ref_grad = central_diff_polish_collaboration(x0, robot.DH, robot.base, robot.cap, curPC);

        s = distMat - ref_grad*theta0;
        l = -ref_grad;
        Sstack = Vcat(Sstack, s);
        Lstack = Vcat(Lstack, l);
    }

    MatrixXd Jac, Diff_Jac, Diff_cfun, Diff;
    Jac = Jacobi(theta0, robot.DH, robot.nlink, c0);
    Diff_Jac = Jac.block(0,0,3,Jac.cols());
    Diff_cfun = Diff_Jac_num_grad(x0, c1, robot.base, robot.cap, robot.tool, PC);
    Diff_cfun = Diff_cfun.block(0, 1, matrix.rows(), 2);
    Diff = Hcat(Diff_Jac,Diff_cfun);

    MatrixXd Aeq, beq;
    Aeq = Diff;
    beq = c1_new - c0 + Diff*x0;



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

bool check_collision_complete_PC_cluster(MatrixXd theta0, Robot robot, MatrixXd PC, MatrixXd PC_idx){
    bool col_flag = false;
    MatrixXd abs_PC_idx = PC_idx.array().abs();
    int max_value = abs_PC_idx.maxCoeff();
    double minDist = INFINITY;
    double dist;
    for (int j=0; j<max_value;j++){
        MatrixXd curPC;
        for (int i = 0; i < PC_idx.rows(); i++) {
            if (PC_idx(i) == j){
                curPC.convervativeResize(curPC.rows() + 1, curPC.cols());
                curPC.row(curPC.rows() - 1) = PC_new.row(i);
            }
        }

        dist = dist_arm_PC(theta0, robot.DH, robot.base, robot.cap, curPC);
        minDist = min(dist, minDist);
        if (minDist <= 0){
            col_flag = true;
            cout << "collision happens" << endl;
            break;
        }
    }
    return col_flag;
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

    if (check_collision_complete_PC_cluster(theta_1,robot,PC_1,PC_idx) == true){
        cout << "P1" << endl;
        valid = false;
    }

    if (checkFeasible(theta_1, PC_1, PC_idx, robot.DH, robot.base, robot.cap) == 0){
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

        if (checkFeasible(theta_step, PC_step, PC_idx, robot.DH, robot.base, robot.cap) == 0) {
            valid = 0;
            cout << "P4" << endl;
            break;
        }

    }
}


bool checkFeasible(MatrixXd theta, MatrixXd PC, MatrixXd PC_idx, MatrixXd DH, MatrixXd base, capsule RoCap[]){
    double PC_minX = PC.col(0).minCoeff();
    double PC_maxX = PC.col(0).maxCoeff();
    double PC_minY = PC.col(1).minCoeff();
    double PC_maxY = PC.col(1).maxCoeff();
    double PC_minZ = PC.col(2).minCoeff();
    double PC_maxZ = PC.col(2).maxCoeff();
    bool feasible = true;

    double offset = 0.2;
    MatrixXd PC_hull = PC.rowwise().any() && (PC_idx.array() < 0).colwise().all();
    DH.col(0) = theta;
    MatrixXd pos, M;
    CapPos(base,DH,RoCap,M,pos);

    for (int i = 0; i < pos.size(); i++) {
        MatrixXd point1 = pos[i].col(0);
        MatrixXd point2 = pos[i].col(1);
        double radius = radius[i];
        if (point1(0) > PC_minX + offset && !inhull(point1.transpose(), PC_hull)) {
            feasible = false;
            return feasible;
        }
        if (point2(0) > PC_minX + offset && !inhull(point2.transpose(), PC_hull)) {
            feasible = false;
            return feasible;
        }
    }

    MatrixXd c_end = ForKine(theta, DH, base, RoCap);
    MatrixXd c_last = pos[7].col(1);
    MatrixXd v = c_end - c_last;
    double normV = v / v.norm();
    MatrixXd normV3d = v / normV;
    MatrixXd c_end1 = c_end - 0.03 * normV3d;

    if (distLinSeg2PC(c_end1, c_last, PC.transpose()) < 0.03 || distLinSeg2PC(c_end, c_last, PC.transpose()) < 0) {
        feasible = false;
        return feasible;
    }

    return feasible;
}


#endif