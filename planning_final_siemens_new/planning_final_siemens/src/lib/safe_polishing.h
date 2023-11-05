#ifndef _SAFEPOLISHING_H_
#define _SAFEPOLISHING_H_

/*
* Safe polishing function definition
* polishing libraries 
* Author -
* Date - 
*/
#define OPTIM_ENABLE_EIGEN_WRAPPERS
//#include <optim/optim.hpp>
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
#include <armadillo>
#include <nlopt.h>
#include <chrono>
// #include "ProblemCFS_qp3d_position.h"
// #include "PlanningProblem_3d_position.h"
using namespace std;
using namespace arma;


// struct myfunc_data {
//     MatrixXd theta_tmp;
//     MatrixXd Jac_null;
//     Robot robot;
//     MatrixXd PC_0;
//     double alpha0;

//     myfunc_data(const MatrixXd& theta_tmp_,
//                 const MatrixXd& Jac_null_,
//                 const Robot& robot_,
//                 const MatrixXd& PC_0_,
//                 double alpha0_)
//         : theta_tmp(theta_tmp_), Jac_null(Jac_null_), robot(robot_), PC_0(PC_0_), alpha0(alpha0_) {
//         // Constructor initializes the members with the provided values
//     }
// };

struct myfunc_data_arma {
    mat theta_tmp;
    mat Jac_null;
    mat DH;
    mat base;
    Robot_arma robot_arma;
    mat PC_0;
    double alpha0;

    myfunc_data_arma(const mat& theta_tmp_,
                     const mat& Jac_null_,
                     const mat& DH_,
                     const mat& base_,
                     const Robot_arma& robot_arma_,
                     const mat& PC_0_,
                     double alpha0_)
        : theta_tmp(theta_tmp_), Jac_null(Jac_null_), DH(DH_), base(base_),
          robot_arma(robot_arma_), PC_0(PC_0_), alpha0(alpha0_) {
        // Initialize the array of pointers to capsule_arma objects
        }
};


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

void addLowUpCons_arma(arma::mat& A, arma::mat& b, arma::mat& x, arma::mat lb, arma::mat ub) {
    /* 
    * this is for quadratic programming 
    * add the lower bound and upper bound for variable x 
    * lb: is the lower bound 
    * ub: is the upper bound
    */
    int dim = x.n_rows + 2;

    // Set the lower bound
    arma::mat Alow = -1 * arma::eye<arma::mat>(dim, dim);

    assert(lb.n_cols == 1); // make sure lb is a column vector
    arma::mat blow = -1 * lb;

    // Set the upper bound
    arma::mat Aup = arma::eye<arma::mat>(dim, dim);

    assert(ub.n_cols == 1);
    arma::mat bup = ub;

    // Concatenate the upper bound condition and lower bound condition together with A, b
    // Concatenate with lower bound
    A = join_vert(A, Alow);
    b = join_vert(b, blow);

    // Concatenate with upper bound
    A = join_vert(A, Aup);
    b = join_vert(b, bup);

    // All is done now.
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

arma::mat Jacobi_arma(arma::mat theta, arma::mat DH, int n, arma::mat p) {
    /* the function calculates the Jacobi matrix (J) of
    * a given point p (in the base coordinate) on the link n
    */
    arma::mat z, r_0, J, TCP_T;
    int nlink = theta.n_rows;
    DH.submat(0, 0, nlink - 1, 0) = theta;
    z = arma::zeros<arma::mat>(3, n); // z(:,i) is the z axis of link i in the base coordinate
    r_0 = arma::zeros<arma::mat>(3, n + 1); // r_0(:,i) is the coordinate of the i-th origin
    J = arma::zeros<arma::mat>(nlink, n); // Jacobian
    TCP_T = arma::eye<arma::mat>(4, 4);

    arma::mat alpha, A, D, q;
    alpha = DH.col(3);
    A = DH.col(2);
    D = DH.col(1);
    q = DH.col(0);

    arma::mat tmp;
    for (int i = 0; i < n; ++i) {
        z.col(i) = TCP_T.submat(0, 2, 2, 2);
        r_0.col(i) = TCP_T.submat(0, 3, 2, 3);
        tmp.resize(4, 4);
        tmp << cos(q(i, 0)) << -sin(q(i, 0)) * cos(alpha(i, 0)) << sin(q(i, 0)) * sin(alpha(i, 0)) << A(i, 0) * cos(q(i, 0))
            << endr
            << sin(q(i, 0)) << cos(q(i, 0)) * cos(alpha(i, 0)) << -cos(q(i, 0)) * sin(alpha(i, 0)) << A(i, 0) * sin(q(i, 0))
            << endr
            << 0 << sin(alpha(i, 0)) << cos(alpha(i, 0)) << D(i, 0)
            << endr
            << 0 << 0 << 0 << 1;
        TCP_T = TCP_T * tmp;
    }
    r_0.col(n) = TCP_T.submat(0, 3, 2, 3);
    // TCP_T=TCP_T*... % up to now, forward kinematics, S0 to S6
    arma::mat up;
    for (int i = 0; i < n; ++i) {
        up = cross(r_0.col(i) - p, z.col(i));
        J.col(i) = join_vert(up, z.col(i));
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

mat readXYZFile_arma(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return mat();
    }

    std::string line;
    std::vector<vec> points;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x, y, z;
        if (!(iss >> x >> y >> z)) {
            std::cerr << "Invalid data format in file: " << filename << std::endl;
            return mat();
        }
        points.push_back(vec({x, y, z}));
    }

    file.close();

    if (points.empty()) {
        std::cerr << "No points found in file: " << filename << std::endl;
        return mat();
    }

    mat pointMatrix(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        pointMatrix.row(i) = points[i].t();
    }

    return pointMatrix;
}


void load_PC(MatrixXd& PC, MatrixXd& PC_idx){
    string sample_file = "parameter/sampled_wp.xyz";
    string xyz_file = "parameter/PC_idx.xyz";
    PC = readXYZFile(sample_file);
    load_PC_idx(PC_idx);
}

void load_PC_arma(mat& PC, mat& PC_idx){
    string sample_file = "parameter/sampled_wp.xyz";
    string xyz_file = "parameter/PC_idx.xyz";
    PC = readXYZFile_arma(sample_file);
    load_PC_idx_arma(PC_idx);
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

mat ForKine_arma(mat theta_ini, mat DH, mat base) {
    mat Msix2tool_arma;
    Msix2tool_arma.resize(4, 4);
    Msix2tool_arma << 0.71906213 << -0.49627146 << -0.48648154 << -0.073912 << endr
      << 0.485629   <<  0.85957094 << -0.15906689 << -0.074794 << endr 
      << 0.49710575 << -0.12187057 <<  0.85908872 <<  0.46376 << endr
      << 0          <<  0          <<  0          <<  1 << endr;

    int nlink = DH.n_rows;
    DH.submat(0, 0, nlink - 1, 0) = theta_ini;
    //assert(all(DH.submat(0, 0, nlink - 1, 0) == theta_ini));

    std::vector<mat> M(nlink + 1);
    mat R, T, JTR;
    M[0].resize(4, 4);
    M[0].eye();

    for (int i = 0; i < nlink; ++i) {
        // Rotation
        R.resize(3, 3);
        R << cos(DH(i, 0)) << -sin(DH(i, 0)) * cos(DH(i, 3)) << sin(DH(i, 0)) * sin(DH(i, 3)) << endr
          << sin(DH(i, 0)) << cos(DH(i, 0)) * cos(DH(i, 3)) << -cos(DH(i, 0)) * sin(DH(i, 3)) << endr
          <<   0 << sin(DH(i, 3)) << cos(DH(i, 3)) << endr;

        // Transformation
        T.resize(3, 1);
        T << DH(i, 2) * cos(DH(i, 0)) << endr
          <<  DH(i, 2) * sin(DH(i, 0)) << endr
            << DH(i, 1) << endr;

        JTR.resize(4, 4);
        JTR.submat(0, 0, 2, 2) = R;
        JTR.submat(0, 3, 2, 3) = T;
        JTR.submat(3, 0, 3, 3).zeros();
        JTR(3, 3) = 1;

        M[i + 1].resize(4, 4);
        M[i + 1] = M[i] * JTR;

        // Additional transformation matrix
        if (i == nlink - 1) {
            M[i + 1] *= Msix2tool_arma;
        }
    }

    // Get the end-effector position
    mat epos = M[nlink].submat(0, 3, 2, 3) + base;
    return epos;
}

double dist_nullspace(MatrixXd x0, MatrixXd theta0, MatrixXd nullspace, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, double alpha){
    MatrixXd theta;
    // MatrixXd zero_tmp(1,1);
    // zero_tmp(0,0) = 0;

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


// double costFunction(unsigned n, const double *x, double *grad, void *data)
// {

//     myfunc_data *params = reinterpret_cast<myfunc_data *>(data);
//     MatrixXd theta_tmp = params->theta_tmp;
//     MatrixXd Jac_null = params->Jac_null;
//     Robot robot = params->robot;
//     MatrixXd PC_0 = params->PC_0;
//     double alpha0 = params->alpha0;

//     MatrixXd xx(3,1);
//     xx << x[0],x[1],x[2];
//     return dist_nullspace(xx, theta_tmp, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
// }



double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
    return (sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])-1);
 }

double dist_nullspace_arma(vec x0, mat theta0, mat nullspace, mat DH, mat base, capsule_arma cap[], tool_arma tl, mat PC, double alpha) {
    mat theta;
    theta = theta0 + alpha * (nullspace * x0);
    int nlink = DH.n_rows;
    DH.submat(0, 0, nlink - 1, 0) = theta;
    double d = INFINITY;
    lineseg_arma* pos = CapPos_real_arma(base, DH, cap, tl);

    mat point1, point2;
    double dist, radius;
    std::vector<double> distList;

    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC_arma(point1, point2, PC.t());
        dist -= radius;
        distList.push_back(dist);
    }

    mat v = ForKine_arma(theta, DH, base) - pos[7].p2;
    double normV = norm(v.col(0));
    //mat normV3d = v / norm;
    mat c_end1 = ForKine_arma(theta, DH, base) - 0.03 * (v / normV);

    double dist_end = distLinSeg2PC_arma(c_end1, c_end1, PC.t()) - 0.03;
    distList.push_back(dist_end);

    // for (int i = 0; i<distList.size();i++)
    // {
    //     cout << distList[i] << endl;
    // }
    auto startIterator = distList.begin() + 3; // Points to the 4th element
    auto endIterator = distList.end();         // Points to the end

    d = *std::min_element(startIterator, endIterator);
    return -d;
}

double costFunction_arma(unsigned n, const double *x, double *grad, void *data)
{

    myfunc_data_arma *params = reinterpret_cast<myfunc_data_arma *>(data);
    mat theta_tmp = params->theta_tmp;
    mat Jac_null = params->Jac_null;
    mat DH = params->DH;
    mat base = params->base;
    Robot_arma robot_arma = params->robot_arma;
    mat PC_0 = params->PC_0;
    double alpha0 = params->alpha0;

    mat xx(3,1);
    xx << x[0] << endr << x[1] << endr << x[2] << endr;
    return dist_nullspace_arma(xx, theta_tmp, Jac_null, DH, base, robot_arma.cap, robot_arma.tl, PC_0, alpha0);
}


double obj_func(const Eigen::VectorXd& vals_inp, Eigen::VectorXd* grad_out, void* opt_data, MatrixXd theta0, MatrixXd nullspace, Robot& robot, MatrixXd PC, double alpha)
{

    const double x = vals_inp(0);
    const double y = vals_inp(1);
    const double z = vals_inp(2);
    MatrixXd input(3,1);
    input << x,y,z;

    const double obj_val = dist_nullspace(input, theta0, nullspace, robot.DH, robot.base, robot.cap, robot.tl, PC, alpha);
            
    return obj_val;
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

void setRobotGoal_arma(Robot_arma& robot, int& nstep1, int& nstep2, int& nwait, int& tfinal, mat& theta_ini, mat& weldTraj, mat& center_point) {
    // Ensure weldTraj dimension is correct
    assert(weldTraj.n_rows == 3);

    // Get pos initial
    mat posini = ForKine_arma(theta_ini, robot.DH, robot.base);
    mat posexe = weldTraj.submat(0, 0, 2, 0);

    // Set robot end effector reference trajectory
    mat diff;
    mat diff1 = center_point - posini;
    mat diff2 = posexe - center_point;
    mat tmp;

    for (int t = 1; t <= nstep1; ++t) {
        tmp = posini + double(t) / nstep1 * diff1;
        robot.goal = join_rows(robot.goal, tmp);
    }

    for (int t = 1; t <= nstep2; ++t) {
        tmp = center_point + double(t) / nstep2 * diff2;
        robot.goal = join_rows(robot.goal, tmp);
    }

    // Wait until the velocity and acceleration approach 0
    for (int t = 1; t <= nwait; ++t) {
        robot.goal = join_rows(robot.goal, posexe);
    }

    // Execute polishing (end-effector located precisely on the weld points)
    int ministep = 1; // Interpolation (augment weld points)
    for (int t = 0; t < tfinal - 1; ++t) {
        diff = weldTraj.submat(0, t + 1, 2, t+1) - weldTraj.submat(0, t, 2, t);
        // Interpolation
        for (int st = 1; st <= ministep; ++st) {
            tmp = weldTraj.submat(0, t, 2, t) + double(st) / ministep * diff;
            robot.goal = join_rows(robot.goal, tmp);
        }
    }
    tfinal = (tfinal - 1) * ministep; // Update tfinal steps
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

bool check_ineq_need_arma(mat theta, mat DH, mat base, capsule_arma cap[], double consider_line){
    int nlink = DH.n_rows;
    DH.col(0) = theta;
    double d = INFINITY;
    lineseg_arma pos[DH.n_rows];
    mat M[DH.n_rows+1];
    CapPos_arma(base, DH, cap, M, pos);
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

        dist = dist_arm_PC_WP(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
        
        distMat.resize(1,1);
        distMat << dist;
        ref_grad = central_diff_polish_collaboration(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);

        s = distMat - ref_grad*x0;
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
    
    QPxset(x, theta0);// set initial value of u;

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

tuple<mat, mat>safetrack_auto_collaboration_arma(mat theta0, mat wp_pos0, Robot_arma robot, mat c1, mat PC, mat PC_idx){
    mat PC_new, wp_pos_new, c1_new, c0, x0;
    c1_new = next_point_WP_arma(wp_pos0, c1, PC);
    double consider_line = 1.2;
    c0 = ForKine_arma(theta0, robot.DH, robot.base);
    x0 = join_cols(theta0,wp_pos0.row(1));
    x0 = join_cols(x0,wp_pos0.row(2));
    assert(robot.nlink == 6);
    double penalty[8] = {5,5,2,1,5,10,50,40};
    double dist;

    mat H(8,8, fill::zeros);
    for (int jnt=0; jnt<8; ++jnt){
        H(jnt,jnt) = penalty[jnt];
        }

    mat f = -H.t()*x0;
    mat Sstack, Lstack, ref_grad, l, s, distMat;
    double max_value = arma::max(arma::abs(PC_idx)).max();

    bool need_flag = check_ineq_need_arma(theta0, robot.DH, robot.base, robot.cap, consider_line);
    if(need_flag){
        for (int j = 1; j <= max_value; j++){
            mat curPC = PC.rows(find(PC_idx == j));
            dist = dist_arm_PC_WP_arma(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
            distMat.resize(1,1);
            distMat << dist;
            ref_grad = central_diff_polish_collaboration_arma(x0, robot.DH, robot.base, robot.cap, robot.tl, curPC);

            s = distMat - ref_grad*x0;
            l = -ref_grad;
            Sstack = join_cols(Sstack, s);
            Lstack = join_cols(Lstack, l);

        }
    }
    mat Jac, Diff_Jac, Diff_cfun, Diff;
    Jac = Jacobi_arma(theta0, robot.DH, robot.nlink, c0);
    Diff_Jac = Jac.submat(0, 0, 2, Jac.n_cols - 1);
    Diff_cfun = Diff_Jac_num_grad_arma(c1, PC, wp_pos0);
    Diff_cfun = Diff_cfun.cols(1, 2);
    Diff = join_rows(Diff_Jac,Diff_cfun);

    mat Aeq, beq;
    Aeq = Diff;
    beq = c1_new - c0 + Diff*x0;
    mat lb = arma::join_vert(robot.lb, arma::vec{-1.0, -1.0});
    mat ub = arma::join_vert(robot.ub, arma::vec{1.0, 1.0});

    addLowUpCons_arma(Lstack, Sstack, theta0, lb, ub);
    quadprogpp::Matrix<double> G, CE, CI;
    quadprogpp::Vector<double> g0, ce0, ci0, x;

    setObjValue_arma(G, g0, H, f); // G is quadratic term and f is the linear term 
    if (Lstack.n_rows==0){// when inequality is not applicable 
        CI.resize(robot.nlink,0);
        ci0.resize(0);
    }

    else{
        mat nLT = -1*Lstack.t(); // negative transpose
        setConstraint_arma(CI, ci0, nLT, Sstack); // set constriants
    }

        // set linear equality cosntraints
    if (beq.n_rows==0){// when equality is not applicable 
        CE.resize(robot.nlink,0);
        ce0.resize(0);
    }
    else{
        mat nAeqT = -1*Aeq.t(); // negative transpose
        setConstraint_arma(CE, ce0, nAeqT, beq); // set cosntriants
    }

    QPxset_arma(x, theta0);  
    double tmp_cost = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    cout << "the temporal cost is " << tmp_cost << endl;
    cout << "valid quadprog is solved" << endl;
    mat x_new, theta_new;
    quad2matrix_arma(x, x_new); // push the solved x value to theta_new
    theta_new = x_new.rows(0, 5);
    wp_pos_new = join_vert(arma::mat(1, 1, arma::fill::zeros),x_new.rows(6,7));
    
    return std::make_tuple(theta_new, wp_pos_new);

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

bool check_collision_complete_PC_cluster_arma(mat theta0, Robot_arma robot, mat PC, mat PC_idx){
    bool col_flag = false;
    double max_value = arma::max(arma::abs(PC_idx)).max();
    double minDist = INFINITY;
    double dist;
    for (int j=1; j<=max_value;j++){
        mat curPC = PC.rows(find(PC_idx == j));
        dist = dist_arm_PC_arma(theta0, robot.DH, robot.base, robot.cap, robot.tl, curPC);
        minDist = min(dist, minDist);
        if (minDist <= -0.0005){
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
    DH.col(0) = theta;
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

bool checkFeasible_arma(mat theta, mat PC, mat PC_idx, mat DH, mat base, capsule_arma RoCap[], tool_arma tl){
    double PC_minX = (PC.col(0)).min();
    bool feasible = true;
    double offset = 0.2;
    mat PC_hull = PC.rows(arma::find(PC_idx < 0));

    DH.col(0) = theta;
    lineseg_arma* pos = CapPos_real_arma(base, DH, RoCap, tl);
    mat point1, point2;
    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        double radius = pos[i].r;

        if (point1(0,0) > (PC_minX + offset) && inhull_arma(point1.t(), PC_hull) == 0) {
        //if (inhull_arma(point1.t(), PC_hull) == 0 && point1(0,0) > (PC_minX + offset)) {
            feasible = false;
            return feasible;
        }
        if (point2(0,0) > (PC_minX + offset) && inhull_arma(point2.t(), PC_hull) == 0) {
            feasible = false;
            return feasible;
        }
    }

    mat c_end = ForKine_arma(theta, DH, base);
    mat c_last = pos[7].p2;
    mat v = c_end - c_last;
    double normV = norm(v.col(0)); 
    mat c_end1 = c_end - 0.03 * (v/normV);
    if (distLinSeg2PC_arma(c_end1, c_last, PC.t()) < 0.03 || distLinSeg2PC_arma(c_end, c_last, PC.t()) < 0) {
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
        return valid;
    }
    if (check_collision_complete_PC_cluster(theta_1,robot,PC_1,PC_idx) == true){
        cout << "P1" << endl;
        valid = false;
        return valid;
    }

    if (checkFeasible(theta_1, PC_1, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0){
        cout << "P2" << endl;
        valid = false;
        return valid;
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

bool checkValid_arma(mat theta_0, mat theta_1, mat wp_pos_0, mat wp_pos_1, mat c1, Robot_arma robot, mat PC_origin, mat PC_idx){
    bool valid = true;
    mat c0 = ForKine_arma(theta_1, robot.DH, robot.base);
    mat PC_0, PC_1;
    tie(PC_0, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_0);
    tie(PC_1, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_1);

    mat diff = c1 - c0;
    if (norm(diff,2) > 0.001){
        cout << "P0" << endl;
        valid = false;
        return valid;
    }
    if (check_collision_complete_PC_cluster_arma(theta_1,robot,PC_1,PC_idx) == true){
        cout << "P1" << endl;
        valid = false;
        return valid;
    }

    if (checkFeasible_arma(theta_1, PC_1, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0){
        cout << "P2" << endl;
        valid = false;
        return valid;
    }

    mat diff_wp_pos = wp_pos_1 - wp_pos_0;
    mat diff_theta = theta_1 - theta_0;
    int steps = 20;
    double minDist = INFINITY;

    mat PC_step = PC_0;
    for (int s = 0; s <= steps; s++){

        mat theta_step = theta_0 + s * diff_theta / steps;
        if (norm(diff_wp_pos,2) > 0) {
            mat wp_pos_step = wp_pos_0 + s * diff_wp_pos/steps;
            tie(PC_step, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_step);
        }

        bool col_flag = check_collision_complete_PC_cluster_arma(theta_step, robot, PC_step, PC_idx);
        if (col_flag == 1) {
            valid = 0;
            cout << "P3" << endl;
            break;
        }

        if (checkFeasible_arma(theta_step, PC_step, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0) {
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
        MatrixXd Jac_new = Jac.topRows(3);
        MatrixXd dtheta = Jac_new.completeOrthogonalDecomposition().pseudoInverse() * dc;
        theta0 = theta0 + dtheta;
    }
    return theta0;
}

mat iterTrack_arma(mat theta0, mat c1, Robot_arma robot){
    for (int i = 0; i < 100; i++){
        mat c0 = ForKine_arma(theta0, robot.DH, robot.base);
        mat dc = c1-c0;
        if(norm(dc,2) < 0.001){
            break;
        }

        mat Jac = Jacobi_arma(theta0, robot.DH, robot.nlink, c0-robot.base);
        mat Jac_new = Jac.rows(0,2);
        mat dtheta = pinv(Jac_new) * dc;
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



void proc_nullspace_opt(MatrixXd theta_cur, MatrixXd wp_pos_cur, Robot robot, MatrixXd PC_origin, MatrixXd PC_idx, MatrixXd c_next, int safe_planner, string nullspace_mode, MatrixXd& theta_pre, MatrixXd& wp_pos_pre, int& success){
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
    double dist = dist_nullspace(x0, theta_tmp, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
    vector<double> d_list;
    d_list.push_back(dist);
    double d_min = dist;
    double alpha = alpha0;
    MatrixXd theta_opt;

    for (int i=1; i<=10; i++) {
        MatrixXd theta_new;
        MatrixXd x_fmin(3,1);
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
            //-----------------------------------------------------------
            // nlopt_opt opt;
            // opt = nlopt_create(NLOPT_LN_COBYLA, 3);
            // // myfunc_data params;
            // // params.theta_tmp = theta_tmp;
            // // params.Jac_null = Jac_null;
            // // params.robot = robot;
            // // params.PC_0 = PC_0;
            // // params.alpha0 = alpha0;
            // myfunc_data params(theta_tmp, Jac_null, robot, PC_0, alpha0);
            // nlopt_set_min_objective(opt, costFunction, &params);
            // nlopt_set_param(opt, "inner_maxeval", 2);
            // nlopt_set_xtol_rel(opt, 1e-4);
            // nlopt_add_equality_constraint(opt, myconstraint, NULL, 1e-8);
            // double x[3] = { 0, 0, 0 };  /* some initial guess */
            // double minf; /* the minimum objective value, upon return */

            // if (nlopt_optimize(opt, x, &minf) < 0) {
            //     cout << "bad" << endl;
            // }
            // else {
            //     cout << x[0] << endl;
            //     cout << x[1] << endl;
            //     cout << x[2] << endl;
            // }
            //-----------------------------------------------------------------
            mat theta_tmp_arma;
            mat Jac_null_arma(Jac_null.rows(),Jac_null.cols());
            mat DH_arma;
            mat base_arma;
            mat PC_0_arma(PC_0.rows(),PC_0.cols());

            theta_tmp_arma.set_size(theta_tmp.rows(), theta_tmp.cols());
            for (int i = 0; i < theta_tmp.rows(); ++i) {
                theta_tmp_arma(i, 0) = theta_tmp(i, 0);
            }

            for (int i = 0; i < Jac_null.rows(); ++i) {
                for (int j = 0; j < Jac_null.cols(); ++j) {
                    Jac_null_arma(i, j) = Jac_null(i, j);
                }
            }

            // if (i == 1){
            //     Jac_null_arma   << -0.1859 <<  -0.0058 <<  -0.0452 << endr
            //                   <<  0.0002  <<  0.2619 <<  -0.0240  << endr
            //                   <<  0.0366 <<  -0.5147 <<   0.0562  << endr
            //                   <<  0.9809  << -0.0174 <<  -0.0031  << endr
            //                   <<  0.0427 <<   0.8161 <<   0.0274  << endr
            //                << -0.0086 <<   0.0126 <<   0.9967 << endr;
            //  }

            // if (i == 2){
            //     Jac_null_arma   << -0.1596 <<  0.0038  << -0.0471 << endr
            //                    << -0.0016  <<  0.2392 <<  -0.0208  << endr
            //                      <<   0.0334 <<  -0.5024 <<   0.0524 << endr
            //                       <<  0.9860 <<  -0.0118 <<  -0.0031 << endr
            //                        << 0.0356 <<   0.8307 <<   0.0248 << endr
            //                     << -0.0071 <<   0.0109 <<  0.9970 << endr;
            // } 

            for (int i = 0; i < PC_0.rows(); ++i) {
                for (int j = 0; j < PC_0.cols(); ++j) {
                    PC_0_arma(i, j) = PC_0(i, j);
                }
            }
            //Jac_null_arma(Jac_null.data(), Jac_null.rows(), Jac_null.cols(), false);
            // PC_0_arma(PC_0.data(), PC_0.rows(), PC_0.cols(), false);
            loadDHbase_arma(DH_arma, base_arma);
            // capsule_arma caps_arma[5];
            // loadCapsules(caps_arma);
            // tool_arma tl_arma;
            // loadTool_arma(tl_arma);
            mat x0_arma(3, 1, fill::zeros);
            vec new_x0 = {0,0,0};
            Robot_arma robot_arma("GP50");
            double dist_arma = dist_nullspace_arma(new_x0, theta_tmp_arma, Jac_null_arma, DH_arma, base_arma, robot_arma.cap, robot_arma.tl, PC_0_arma, alpha0);
            nlopt_opt opt;
            opt = nlopt_create(NLOPT_LN_COBYLA, 3);

            // opt = nlopt_create(NLOPT_GN_AGS, 3);
            myfunc_data_arma params(theta_tmp_arma, Jac_null_arma, DH_arma, base_arma, robot_arma, PC_0_arma, alpha0);
            nlopt_set_min_objective(opt, costFunction_arma, &params);
            //nlopt_set_param(opt, "inner_maxeval", 5);
            nlopt_set_xtol_rel(opt, 1e-3);

            double lb[3] = { -1.0, -1.0, -1.0 }; //lower bounds
            double ub[3] = {1.0 , 1.0, 1.0}; //upper bounds
            nlopt_set_lower_bounds(opt,lb);
            nlopt_set_upper_bounds(opt,ub);
            nlopt_set_ftol_rel(opt, 1e-3);
            nlopt_set_maxeval(opt, 10);
            nlopt_add_equality_constraint(opt, myconstraint, NULL, 1e-3);
            double x[3] = { 0, 0, 0 };  /* some initial guess */
            double minf; /* the minimum objective value, upon return */
            // //double dist_arma2 = dist_nullspace_arma(x, theta_tmp_arma, Jac_null_arma, DH_arma, base_arma, caps_arma, tl_arma, PC_0_arma, alpha0);
            // //cout << "dist arma" << dist_arma << endl;
            if (nlopt_optimize(opt, x, &minf) < 0) {
                cout << "bad" << endl;
            }
            else {
                cout << x[0] << endl;
                cout << x[1] << endl;
                cout << x[2] << endl;
                x_fmin << x[0],x[1],x[2];
                cout << nlopt_get_numevals(opt) << endl;
                cout << nlopt_get_xtol_rel(opt) << endl;
                cout << nlopt_get_ftol_rel(opt) << endl;


            }
        }
            theta_new = theta_tmp + alpha*(Jac_null*x_fmin);
            theta_new = iterTrack(theta_new,c1,robot);
            double dist_new = dist_nullspace(x0, theta_new, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
            if (dist_new >= d_list.back()){
                x0.setZero();
            }

            theta_tmp = theta_new;
            if(dist_new < d_min && checkValid(theta_cur, theta_new, wp_pos_cur, wp_pos_cur, c1, robot, PC_origin, PC_idx) == 1){
                d_min = dist_new;
                theta_opt = theta_new;
            }

            d_list.push_back(dist_new);
            if(i>5 && d_min < 0){
                break;
            }

            // obj_fun f = [theta_tmp_arma, Jac_null_arma, DH_arma, base_arma, &caps_arma, &tl_arma, PC_0_arma, alpha0](vec &x) -> double
        	// {
		    //     return dist_nullspace_arma(x, theta_tmp_arma, Jac_null_arma, DH_arma, base_arma, caps_arma, tl_arma, PC_0_arma, alpha0);
	        // };

            // // // vec new_x0 = {1,0,0};

            // // // double test_result = f(new_x0);
            // // // cout <<"new results"<< test_result<<endl;
            // vec lb = {-1,-1,-1};
            // vec ub = {1,1,1};
            // auto c = [](vec &x) -> vec
	        // {
            //     return {norm(x)-1, -norm(x)+1};
	        // };
        	// options opt;
	        // opt.algo = Powell;
            // //opt.tolerance = 1e-6;
            // opt.max_ite = 10;
            // opt.tolerance = 0.01;
	        // auto result = sci_arma::fmincon(f, new_x0, lb, ub, c, opt);
            // cout << result << endl;
            // -----------------------------------------------------------------------
        }
        if (theta_opt.rows() > 0){
            theta_pre = theta_opt;
            wp_pos_pre = wp_pos_cur;
            success = 1;
        }

                

}

void proc_nullspace_opt_arma(mat theta_cur, mat wp_pos_cur, Robot_arma robot, mat PC_origin, mat PC_idx, mat c_next, int safe_planner, string nullspace_mode, mat& theta_pre, mat& wp_pos_pre, int& success){
    double alpha0 = 0.1;
    mat theta_tmp = theta_cur;
    mat wp_pos_tmp = wp_pos_cur;
    mat x0(3,1,fill::zeros);
    mat c0 = ForKine_arma(theta_tmp, robot.DH, robot.base);
    mat PC_0;
    tie(PC_0, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_cur);
    mat c1 = next_point_WP_arma(wp_pos_cur, c_next, PC_origin);
    theta_tmp = iterTrack_arma(theta_cur, c1, robot);
    mat Jac = Jacobi_arma(theta_tmp, robot.DH, robot.nlink, c0 - robot.base);
    mat Jac_new = Jac.rows(0,2);

    mat Jac_null = null(Jac_new);

    double dist = dist_nullspace_arma(x0, theta_tmp, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);

    vector<double> d_list;
    d_list.push_back(dist);
    double d_min = dist;
    double alpha = alpha0;
    mat theta_opt;

    for (int i=1; i<=20; i++) {
        mat theta_new;
        mat x_fmin(3,1);
        if(dist < 0 && safe_planner == 0 && checkValid_arma(theta_cur, theta_tmp, wp_pos_cur, wp_pos_cur, c1, robot, PC_origin, PC_idx) == 1){
            theta_opt = theta_tmp;
            break;
        }
        c0 = ForKine_arma(theta_tmp, robot.DH, robot.base);
        Jac = Jacobi_arma(theta_tmp, robot.DH, robot.nlink, c0-robot.base);
        Jac_new = Jac.rows(0,2);
        Jac_null = null(Jac_new);

        if (nullspace_mode == "cmaes"){
            continue;
        }
        else{
            //-----------------------------------------------------------------
            mat x0_arma(3, 1, fill::zeros);
            vec new_x0 = {0,0,0};
            double dist_arma = dist_nullspace_arma(new_x0, theta_tmp, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
            nlopt_opt opt;
            opt = nlopt_create(NLOPT_LN_COBYLA, 3);
            // opt = nlopt_create(NLOPT_GN_AGS, 3);
            myfunc_data_arma params(theta_tmp, Jac_null, robot.DH, robot.base, robot, PC_0, alpha0);
            nlopt_set_min_objective(opt, costFunction_arma, &params);
            //nlopt_set_param(opt, "inner_maxeval", 5);
            //nlopt_set_xtol_rel(opt, 5e-3);
            //nlopt_set_xtol_rel(opt, 5e-2);
            double lb[3] = { -1.0, -1.0, -1.0 }; //lower bounds
            double ub[3] = {1.0 , 1.0, 1.0}; //upper bounds
            nlopt_set_lower_bounds(opt,lb);
            nlopt_set_upper_bounds(opt,ub);
            //nlopt_set_ftol_rel(opt, 1e-3);
            //nlopt_set_ftol_rel(opt, 5e-2);
            nlopt_set_maxeval(opt, 35);
            //nlopt_add_equality_constraint(opt, myconstraint, NULL, 1e-3);
            //nlopt_add_equality_constraint(opt, myconstraint, NULL, 1e-2);

            double x[3] = { 0, 0, 0 };  /* some initial guess */
            double minf; /* the minimum objective value, upon return */
            if (nlopt_optimize(opt, x, &minf) < 0) {
                cout << "bad" << endl;
            }
            else {
                // cout << x[0] << endl;
                // cout << x[1] << endl;
                // cout << x[2] << endl;
                x_fmin << x[0] << endr
                        << x[1] << endr
                        << x[2] << endr;
                // cout << nlopt_get_numevals(opt) << endl;
                // cout << nlopt_get_xtol_rel(opt) << endl;
                // cout << nlopt_get_ftol_rel(opt) << endl;
            }
            //x_fmin = fmincon_x.col(ft);
            ft++;
        }
            theta_new = theta_tmp + alpha*(Jac_null*x_fmin);
            //cout << "theta_new" << theta_new << endl;

            theta_new = iterTrack_arma(theta_new,c1,robot);
            //cout << "theta_new" << theta_new << endl;

            double dist_new = dist_nullspace_arma(x0, theta_new, Jac_null, robot.DH, robot.base, robot.cap, robot.tl, PC_0, alpha0);
            //cout << "dist_new" << dist_new << endl;

            if (dist_new >= d_list.back()){
                x0.zeros();
            }

            theta_tmp = theta_new;
            if(dist_new < d_min && checkValid_arma(theta_cur, theta_new, wp_pos_cur, wp_pos_cur, c1, robot, PC_origin, PC_idx) == 1){
                d_min = dist_new;
                theta_opt = theta_new;
            }

            d_list.push_back(dist_new);
            if(i>5 && d_min < 0){
                break;
            }

        }
        if (theta_opt.n_rows > 0){
            theta_pre = theta_opt;
            wp_pos_pre = wp_pos_cur;
            success = 1;
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


cout << (tfinal + nstep1 + nstep2 + nwait) << endl;
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
                proc_nullspace_opt(theta_cur, wp_pos_cur, robot, PC_origin, PC_idx, c_next, safe_planner, nullspace_mode, theta_pre, wp_pos_pre, success);
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


void safe_polish_procedure_arma(mat& safe_theta, mat& safe_traj, mat& safe_wp_pos, Robot_arma robot, int nstep1, int nstep2, int nwait, int tfinal, 
                mat &theta_pre, mat &c_pre, double thres, mat &M_PC_0, mat &wp_pos_pre, string track_mode_main, 
                string track_solver, mat PC_origin, mat PC_idx, string explore_mode, string nullspace_mode){
mat c_next, diff;
bool col_flag;
mat theta_new, wp_pos_new;
int repeat_cnt;

vector<double> curLog(5);
vector<vector<double>> log;

for (int t=0; t<tfinal + nstep1 + nstep2 + nwait; ++t){
    cout << "---------- time step " << t << " ----------" << endl;
    c_next = robot.goal.submat(0, t, 2, t);
    mat inv = M_PC_0.i();
    c_next = setVertice_arma(c_next.t(), inv).t();
    curLog[0] = t;
    col_flag = true;

    mat theta_cur = theta_pre;
    mat c_cur = c_pre;
    mat wp_pos_cur = wp_pos_pre;
    mat PC_new;

    int iteration_cnt = 1;
    while(iteration_cnt <= 20){
        iteration_cnt++;
        if (track_mode_main == "robot")
        {
            if (track_solver == "ICOP")
            {
                cout << "****************Robot ICOP Tracking*****************" << endl;
                continue;
                 // tie(theta_new, wp_pos_new) = safetrack_auto_robot(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
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
            tie(theta_new, wp_pos_new) = safetrack_auto_collaboration_arma(theta_pre, wp_pos_pre, robot, c_next, PC_origin, PC_idx);
        }       

        if (theta_new.size() == 0)
            continue;

        mat c_new = ForKine_arma(theta_new, robot.DH, robot.base);
        tie(PC_new, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_new);
        mat c_target = next_point_WP_arma(wp_pos_new, c_next, PC_origin);
        c_pre = c_new;
        theta_pre = theta_new;
        wp_pos_pre = wp_pos_new;
        
        col_flag = check_collision_complete_PC_cluster_arma(theta_new, robot, PC_new, PC_idx);
        diff = c_target - c_pre;
        if (norm(diff,2) <= thres && col_flag == 0){
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

    mat c1 = next_point_WP_arma(wp_pos_pre, c_next, PC_origin);
    if (checkValid_arma(theta_cur, theta_pre, wp_pos_cur, wp_pos_pre, c1, robot, PC_origin, PC_idx) == 0){
        safe_planner = 1;
    }

    if (safe_planner || iteration_cnt > 20 || checkFeasible_arma(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap, robot.tl) == 0){
        if (explore_mode == "None") {
            stop_flag = 1;
        }
        else {
            need_sample = 1;
            success = 0;

            if (explore_mode == "NSopt") {
                proc_nullspace_opt_arma(theta_cur, wp_pos_cur, robot, PC_origin, PC_idx, c_next, safe_planner, nullspace_mode, theta_pre, wp_pos_pre, success);
            }
            
            if (success == 0){
                stop_flag =1;
            }
            else {
                tie(PC_new, ignore, ignore, ignore) = processPC_arma(PC_origin, wp_pos_pre);
                if (checkFeasible_arma(theta_pre, PC_new, PC_idx, robot.DH, robot.base, robot.cap,robot.tl) == 0) {
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
    safe_traj = join_rows(safe_traj, c_pre);
    safe_theta = join_rows(safe_theta, theta_pre);
    safe_wp_pos = join_rows(safe_wp_pos, wp_pos_pre);

    cout << "safe_traj" << safe_traj <<endl;
    cout << "safe_theta" << safe_theta << endl;
    cout << "safe_wp_pos" << safe_wp_pos << endl;
    cout << "theta_pre" << theta_pre << endl;

    }

}


void safe_polish(MatrixXd &start2exe_traj, MatrixXd &execution_traj, MatrixXd &exit_traj){

string name = "GP50";
Robot robot(name);
Robot_arma robot_arma(name);
cout << "DH:\n" << robot.DH << endl;
cout << robot_arma.DH;
cout << "nliink:\n" << robot.nlink << endl;




// load_PC
MatrixXd PC, PC_idx;
load_PC(PC, PC_idx);

mat PC_arma, PC_idx_arma;
load_PC_arma(PC_arma, PC_idx_arma);


MatrixXd arr_axis3;
loadWeldTraj(arr_axis3);
mat arr_axis3_arma;
loadWeldTraj_arma(arr_axis3_arma);

mat arr_axis3T_arma = arr_axis3_arma.t();
MatrixXd arr_axis3T = arr_axis3.transpose();


MatrixXd abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top;

loadWorkpieceSetting(abc, planePoints, point_anchor_axis3, weld_bottom, weld_in, weld_left, weld_out, weld_right, weld_top);
mat abc_arma, planePoints_arma, point_anchor_axis3_arma, weld_bottom_arma, weld_in_arma, weld_left_arma, weld_out_arma, weld_right_arma, weld_top_arma;
loadWorkpieceSetting_arma(abc_arma, planePoints_arma, point_anchor_axis3_arma, weld_bottom_arma, weld_in_arma, weld_left_arma, weld_out_arma, weld_right_arma, weld_top_arma);



MatrixXd weld_in_T = weld_in.transpose();
mat weld_in_T_arma = weld_in_arma.t();

MatrixXd in_center = calculateRowMeans(weld_in_T);
mat in_center_arma = mean(weld_in_T_arma,1);

MatrixXd radius = weld_in_T.colwise() - in_center.col(0);
mat radius_arma = weld_in_T_arma.each_col() - in_center_arma;

VectorXd vecnorm = radius.colwise().norm();
radius = 0.035*radius;
radius.transposeInPlace();
MatrixXd temp = radius.array().colwise()/vecnorm.array();
weld_in = weld_in - temp;

mat vecnorm_arma = arma::sqrt(arma::sum(arma::square(radius_arma), 0));
radius_arma = 0.035*radius_arma;
mat temp_arma = radius_arma.each_row()/vecnorm_arma.as_row();
weld_in_arma = weld_in_arma - temp_arma.t();

//weld_in_arma = weld_in_arma - 0.035 * radius_arma.each_row() / arma::sqrt(arma::sum(arma::square(radius_arma), 0));


MatrixXd temp1(1,3);
temp1 << 0,0.02,0;
MatrixXd temp2(1,3);
temp2 << 0,0,0.02;

mat temp1_arma(1,3);
mat temp2_arma(1,3);
temp1_arma = {{0,0.02,0}};
temp2_arma = {{0,0,0.02}};

weld_right = weld_right.rowwise() + temp1.row(0);
weld_left = weld_left.rowwise() - temp1.row(0);
weld_bottom = weld_bottom.rowwise() + temp2.row(0);
weld_top = weld_top.rowwise() - temp2.row(0);
arr_axis3 = weld_right.transpose();

weld_right_arma.each_row() += temp1_arma;
weld_left_arma.each_row() -= temp1_arma;
weld_bottom_arma.each_row() += temp2_arma;
weld_top_arma.each_row() -= temp2_arma;
arr_axis3_arma = weld_right_arma.t();

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

mat theta_init_default_arma(6,1), wp_pos_init_default_arma(3,1);
loadSafePolishingSetting_arma(alphaY_limit, alphaZ_limit, thres, nstep1, nstep2, nwait, theta_init_default_arma, wp_pos_init_default_arma);

vector<Vector3d> wp_pos_list;

double step = 0.1;

// for (double alphaY = -alphaY_limit; alphaY <= alphaY_limit; alphaY += step) {
//     for (double alphaZ = -alphaZ_limit; alphaZ <= alphaZ_limit; alphaZ += step){
//             Vector3d wp_pos_cur(0.0, alphaY + 0.5, alphaZ);
//             mat wp_pos_cur_arma;
//             wp_pos_list.push_back(wp_pos_cur);
//             MatrixXd PC_processed, M_PC_processed;
//             MatrixXd base_point_processed, center_point_processed;
//             tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_cur);
//             MatrixXd arr_cur = setVertice(arr_axis3.transpose(), M_PC_processed).transpose();
//     }

// }

double joint2_limit = 0.2;
double joint3_limit = 0.2;
double joint4_limit = 0.2;
//vector<VectorXd> theta_init_list;

// for (double joint2 = -joint2_limit; joint2 <= joint2_limit; joint2 += step) {
//     for (double joint3 = -joint3_limit; joint3 <= joint3_limit; joint3 += step) {
//         for (double joint4 = -joint4_limit; joint4 <= joint4_limit; joint4 += step) {
//             VectorXd theta_init(6);
//             theta_init << 0.0, -2.3 + joint2, 0.5 + joint3, 0.0 + joint4, 0.5, -2.0;
//             theta_init_list.push_back(theta_init);
//             MatrixXd PC_processed, M_PC_processed;
//             VectorXd base_point_processed, center_point_processed;
//             tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_init_default);
//             MatrixXd arr_cur = setVertice(arr_axis3.transpose(), M_PC_processed);
//         }

//     }

// }

MatrixXd theta_init = theta_init_default;
MatrixXd wp_pos_init = wp_pos_init_default;
MatrixXd PC_processed, M_PC_processed;
MatrixXd base_point_processed, center_point_processed;
tie(PC_processed, M_PC_processed, base_point_processed, center_point_processed) = processPC(PC, wp_pos_init);
MatrixXd arr = setVertice(arr_axis3.transpose(), M_PC_processed).transpose();

mat theta_init_arma = theta_init_default_arma;
mat wp_pos_init_arma = wp_pos_init_default_arma;
mat PC_processed_arma, M_PC_processed_arma;
mat base_point_processed_arma, center_point_processed_arma;
tie(PC_processed_arma, M_PC_processed_arma, base_point_processed_arma, center_point_processed_arma) = processPC_arma(PC_arma, wp_pos_init_arma);
mat arr_arma = setVertice_arma(arr_axis3_arma.t(), M_PC_processed_arma).t();


int tfinal = arr_arma.n_cols;

setRobotGoal(robot, nstep1, nstep2, nwait, tfinal, theta_init, arr, center_point_processed);
setRobotGoal_arma(robot_arma, nstep1, nstep1, nwait, tfinal, theta_init_arma, arr_arma, center_point_processed_arma);

MatrixXd theta_pre, wp_pos_pre, c_pre, M_PC_0;
theta_pre = theta_init;
wp_pos_pre = wp_pos_init;
c_pre = ForKine(theta_init, robot.DH, robot.base, robot.cap);
M_PC_0 = M_PC_processed;

mat theta_pre_arma, wp_pos_pre_arma, c_pre_arma, M_PC_0_arma;
theta_pre_arma = theta_init_arma;
wp_pos_pre_arma = wp_pos_init_arma;
c_pre_arma = ForKine_arma(theta_init_arma, robot_arma.DH, robot_arma.base);
M_PC_0_arma = M_PC_processed_arma;


MatrixXd safe_traj = c_pre;
MatrixXd safe_theta = theta_pre;
MatrixXd safe_wp_pos = wp_pos_pre;

mat safe_traj_arma = c_pre_arma;
mat safe_theta_arma = theta_pre_arma;
mat safe_wp_pos_arma = wp_pos_pre_arma;

// safe_polish_procedure(safe_theta, safe_traj, safe_wp_pos, robot, 
// nstep1, nstep2, nwait, tfinal, theta_pre, c_pre, thres, M_PC_0, 
// wp_pos_pre, track_mode_main, track_solver, PC, PC_idx, explore_mode, nullspace_mode);

safe_polish_procedure_arma(safe_theta_arma, safe_traj_arma, safe_wp_pos_arma, robot_arma, 
nstep1, nstep2, nwait, tfinal, theta_pre_arma, c_pre_arma, thres, M_PC_0_arma, 
wp_pos_pre_arma, track_mode_main, track_solver, PC_arma, PC_idx_arma, explore_mode, nullspace_mode);



std::ofstream file("safe_theta_test.txt"); // Open a file for writing

if (file.is_open()) {
    for (int i = 0; i < safe_theta_arma.n_rows; ++i) {
        for (int j = 0; j < safe_theta_arma.n_cols; ++j) {
            file << safe_theta_arma(i, j) << " ";
        }
        file << "\n";
    }

    file.close();
    std::cout << "Matrix written to safe_theta_test.txt" << std::endl;
} else {
    std::cerr << "Unable to open file." << std::endl;
}

std::ofstream file1("safe_wp_pos_test.txt"); // Open a file for writing

if (file1.is_open()) {
    for (int i = 0; i < safe_wp_pos_arma.n_rows; ++i) {
        for (int j = 0; j < safe_wp_pos_arma.n_cols; ++j) {
            file1 << safe_wp_pos_arma(i, j) << " ";
        }
        file1 << "\n";
    }
    file1.close();
    std::cout << "Matrix written to safe_wp_pos.txt" << std::endl;
} else {
    std::cerr << "Unable to open file." << std::endl;
}



// MatrixXd safe_theta_T = safe_theta.transpose();
// vector<VectorXd> uniqueRows;

// for (int i = 0; i < safe_theta_T.rows(); ++i) {
//     VectorXd currentRow = safe_theta_T.row(i);

//     bool isUnique = true;
//     for (const auto& row : uniqueRows) {
//       if (currentRow.isApprox(row)) {
//         isUnique = false;
//         break;
//       }
//     }
//     if (isUnique) {
//       uniqueRows.push_back(currentRow);
//     }
// }

// MatrixXd safe_theta_unique(uniqueRows.size(), safe_theta_T.cols());

// for (int i = 0; i < uniqueRows.size(); ++i) {
//     safe_theta_unique.row(i) = uniqueRows[i];
//     }   

// MatrixXd safe_theta_unique_final = safe_theta_unique.transpose();
}











#endif