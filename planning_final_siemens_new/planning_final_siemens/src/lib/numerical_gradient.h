#ifndef _NGrad_
#define _NGrad_

/*
* Numerical gradient library 
* Use 2nd order central difference 
* Author - 
* Date - 
*/

#include "safe_polishing.h"
#include <Eigen/Dense>
#include "structure.h"
#include "distance_constraint_3d.h"
#include <cmath>
#include <armadillo>

using namespace Eigen;
using namespace std;
using namespace arma;


MatrixXd setVertice(MatrixXd v, MatrixXd TransM){
    MatrixXd RotM = TransM.block(0,0,3,3);
    MatrixXd TransVec = TransM.block(0,3,3,1).transpose();
    int n_row = v.rows();
    MatrixXd temp_mat = TransVec.replicate(n_row,1);
    MatrixXd result = v*RotM.transpose() + temp_mat;

    return result;
}

arma::mat setVertice_arma(const arma::mat& v, const arma::mat& TransM) {
    arma::mat RotM = TransM.submat(0, 0, 2, 2);
    arma::mat TransVec = TransM.submat(0, 3, 2, 3).t();
    int n_row = v.n_rows;
    arma::mat temp_mat = repmat(TransVec, n_row, 1);
    arma::mat result = v * RotM.t() + temp_mat;

    return result;
}

tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> processPC(const MatrixXd PC, const MatrixXd PC_pos){
    MatrixXd M1;
    // M1.resize(4,4);
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

    MatrixXd center_point = calculateColMeans(PC_x0);

    MatrixXd mean_point = calculateColMeans(PC).transpose();

    MatrixXd base_point = mean_point.transpose();

    base_point(2) = PC.col(2).minCoeff() - 0.3;

    base_point = setVertice(base_point, M1).transpose();


    MatrixXd M_PC = M1;
    if (PC_pos.size() >= 3) {
        Vector3d elu = M1.block<3, 3>(0, 0).eulerAngles(0, 1, 2);

        for (int i = 0; i < 3; ++i) {
            MatrixXd v(1,3);
            v.setZero();
            v(0,i) = 1;
            if (elu(i) == 0) {
                continue;
            }
            M_PC = mtx_rotate(-1*v*elu(i), base_point) * M_PC;
        }

        for (int i = 0; i < 3; ++i) {
            MatrixXd v(1,3);
            v.setZero();
            v(0,i) = 1;
            if (PC_pos(i,0) == 0) {
                continue;
            }
            M_PC = mtx_rotate(v*PC_pos(i,0), base_point) * M_PC;

        }
    }

    MatrixXd PC_modified = PC;
    MatrixXd center_point_modified = center_point;

    PC_modified = setVertice(PC, M_PC);

    center_point_modified = setVertice(center_point, M_PC).transpose();

    return std::make_tuple(PC_modified, M_PC, base_point, center_point_modified);

}

tuple<mat, mat, mat, mat> processPC_arma(arma::mat PC, const arma::mat PC_pos) {
    // Rotation: ±200°
    // Tilt: ±135°
    // Siemens 
    arma::mat M1 = {{0.91647976, -0.00729021, 0.40001462, 1.40608102},
                    {0.02281488, 0.99915928, -0.03406199, -0.02509452},
                    {-0.39943, 0.0403434, 0.91587558, 1.04437673},
                    {0, 0, 0, 1}};
    MatrixXd M1_eigen;
    loadM1(M1_eigen);
    arma::vec x = PC.col(0);
    arma::uvec idx = arma::find(x == 0);
    arma::mat PC_x0 = PC.rows(idx);
    arma::mat center_point = arma::mean(PC_x0, 0);
    arma::mat mean_point = arma::mean(PC, 0).t();
    arma::mat base_point = mean_point.t();
    base_point(0,2) = arma::min(PC.col(2)) - 0.3;

    base_point = setVertice_arma(base_point, M1).t();

    arma::mat M_PC = M1;
    if (PC_pos.n_elem >= 3) {
        //arma::vec elu;
        Vector3d elu = M1_eigen.block<3, 3>(0, 0).eulerAngles(0, 1, 2);

        for (int i = 0; i < 3; ++i) {
            arma::mat v = arma::zeros<arma::mat>(1,3);
            v(0,i) = 1;
            if (std::abs(elu(i)) < 1e-9) {
                continue;
            }
            M_PC = mtx_rotate_arma(-1*elu(i)*v, base_point) * M_PC;

        }

        for (int i = 0; i < 3; ++i) {
            arma::mat v = arma::zeros<arma::mat>(1,3);
            v(0,i) = 1;
            if (PC_pos(i,0) == 0) {
                continue;
            }
            M_PC = mtx_rotate_arma(v * PC_pos(i,0), base_point) * M_PC;

        }
    }

    PC = setVertice_arma(PC, M_PC);
    center_point = setVertice_arma(center_point, M_PC).t();
    return std::make_tuple(PC, M_PC, base_point, center_point);

}

MatrixXd next_point_WP(MatrixXd wp_pos, MatrixXd c1, MatrixXd PC_origin) {
    MatrixXd curPC;
    MatrixXd M_PC;
    tie(curPC, M_PC, ignore, ignore) = processPC(PC_origin, wp_pos);
    MatrixXd c_next = setVertice(c1.transpose(), M_PC).transpose();
    return c_next;
}

mat next_point_WP_arma(mat wp_pos, mat c1, mat PC_origin) {
    mat curPC;
    mat M_PC;
    tie(curPC, M_PC, ignore, ignore) = processPC_arma(PC_origin, wp_pos);
    mat c_next = setVertice_arma(c1.t(), M_PC).t();
    return c_next;
}



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

double dist_arm_PC_WP(MatrixXd x0, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, int check = 0){
    MatrixXd theta, PC_Processed;
    MatrixXd wp_pos(3,1);
    theta = x0.block(0,0,6,x0.cols());
    wp_pos << 0.0, x0(6,0), x0(7,0);

    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;
    tie(PC_Processed, ignore, ignore, ignore) = processPC(PC, wp_pos);

    lineseg* pos = CapPos_real(base, DH, cap, tl);
    PC_Processed.transposeInPlace();
    MatrixXd point1, point2;
    double dist,radius;
    vector<double> distList(sizeof(pos));
    
    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC(point1, point2, PC_Processed);
        dist -= radius;

        distList[i] = dist;
    }

    d = *std::min_element(distList.begin(),distList.end());
    return d;
}

double dist_arm_PC_WP_arma(mat x0, mat DH, mat base, capsule_arma cap[], tool_arma tl, mat PC, int check = 0){
    mat theta, PC_Processed;
    mat wp_pos(3,1);
    theta = x0.submat(0, 0, 5, x0.n_cols - 1);
    wp_pos << 0.0 << endr 
            << x0(6,0) << endr 
            << x0(7,0) << endr;

    int nlink = DH.n_rows;
    DH.col(0) = theta;
    double d = INFINITY;
    tie(PC_Processed, ignore, ignore, ignore) = processPC_arma(PC, wp_pos);

    lineseg_arma* pos = CapPos_real_arma(base, DH, cap, tl);
    inplace_trans(PC_Processed);
    mat point1, point2;
    double dist,radius;
    vector<double> distList;
    
    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;

        dist = distLinSeg2PC_arma(point1, point2, PC_Processed);
        dist -= radius;

        distList.push_back(dist);
    }

    d = *std::min_element(distList.begin(),distList.end());
    return d;
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

mat central_diff_polish_collaboration_arma(mat x0, mat DH, mat base, capsule_arma cap[], tool_arma tl, mat PC){
    int dim = x0.n_rows;
    double y = dist_arm_PC_WP_arma(x0, DH, base, cap, tl, PC);
    mat grad(1, dim);

    double eps = 1e-5;
    mat x = x0;
    mat xp = x0;
    assert(x0.n_cols == 1);

    double yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = dist_arm_PC_WP_arma(x0, DH, base, cap, tl, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = dist_arm_PC_WP_arma(x0, DH, base, cap, tl, PC);
        grad(0,i) = (yhi - ylo)/eps;
        xp(i,0) = x(i,0);
    }
    return grad; 
}

MatrixXd Diff_Jac_num_grad(MatrixXd c1, MatrixXd PC, MatrixXd x){
    int dim = x.rows();
    MatrixXd y = next_point_WP(x, c1, PC);
    MatrixXd grad(y.rows(), dim);

    double eps = 1e-5;
    MatrixXd xp = x;
    assert(x.cols() == 1);

    MatrixXd yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = next_point_WP(xp, c1, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = next_point_WP(xp, c1, PC);
        grad.col(i) = (yhi - ylo)/eps;
        xp(i,0) = x(i,0);
    }

    return grad; 
}

mat Diff_Jac_num_grad_arma(mat c1, mat PC, mat x){
    int dim = x.n_rows;
    mat y = next_point_WP_arma(x, c1, PC);
    mat grad(y.n_rows, dim);

    double eps = 1e-5;
    mat xp = x;
    assert(x.n_cols == 1);

    mat yhi, ylo;
    for (int i=0; i<dim; ++i){
        xp(i,0) = x(i,0)+eps/2;
        yhi = next_point_WP_arma(xp, c1, PC);
        xp(i,0) = x(i,0)-eps/2;
        ylo = next_point_WP_arma(xp, c1, PC);
        grad.col(i) = (yhi - ylo)/eps;
        xp(i,0) = x(i,0);
    }

    return grad; 
}


#endif