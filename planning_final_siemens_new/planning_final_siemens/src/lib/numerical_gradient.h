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

using namespace Eigen;
using namespace std;

// MatrixXd setVertice(MatrixXd v, MatrixXd TransM) {
//     MatrixXd result = v * TransM.block(0, 0, 3, 3).transpose();
//     result.rowwise() += TransM.block(0, 3, 3, 1).transpose().replicate(v.rows(), 1);
//     return result;
// }

// MatrixXd setVertice(MatrixXd v, MatrixXd TransM) {
//     MatrixXd v_transformed = v*TransM.block(0,0,3,3).transpose();
//     v_transformed.rowwise() += TransM.block(0,3,3,1).transpose();
//     return v_transformed;
// }

// Eigen::MatrixXd setVertice(const Eigen::MatrixXd& v, const Eigen::MatrixXd& TransM) {
//     // Get the dimensions of the input matrix 'v'
//     int num_rows = v.rows();
//     int num_cols = v.cols();

//     // Extract the rotation matrix from 'TransM'
//     Eigen::MatrixXd rotation = TransM.block(0, 0, 3, 3);

//     // Extract the translation vector from 'TransM'
//     Eigen::VectorXd translation = TransM.block(0, 3, 3, 1).transpose();

//     // Initialize the output matrix 'result' with the same dimensions as 'v'
//     Eigen::MatrixXd result(num_rows, num_cols);

//     // Perform the transformation for each row of 'v'
//     for (int i = 0; i < num_rows; ++i) {
//         Eigen::VectorXd row = v.row(i); // Extract each row of 'v' as a vector
//         Eigen::VectorXd transformed_row = row.transpose() * rotation; // Apply rotation
//         transformed_row += translation; // Apply translation

//         // Store the transformed row in the output matrix
//         result.row(i) = transformed_row;
//     }

//     return result;
// }

MatrixXd setVertice(MatrixXd v, MatrixXd TransM){
    MatrixXd RotM = TransM.block(0,0,3,3);
    MatrixXd TransVec = TransM.block(0,3,3,1).transpose();

    // MatrixXd temp_mat = TransVec;
    // cout << "temp_mat" << temp_mat << endl;

    int n_row = v.rows();
    // if(n_row > 1){
    //     for (int i =1; i<n_row; ++i){
    //         temp_mat = Vcat(temp_mat,TransVec);
    //     }
    // }
    MatrixXd temp_mat = TransVec.replicate(n_row,1);
    MatrixXd result = v*RotM.transpose() + temp_mat;

    return result;
}

tuple<MatrixXd, MatrixXd, MatrixXd, MatrixXd> processPC(const MatrixXd& PC, const MatrixXd& PC_pos){
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

MatrixXd next_point_WP(MatrixXd wp_pos, MatrixXd c1, MatrixXd PC_origin) {
    MatrixXd curPC;
    MatrixXd M_PC;
    tie(curPC, M_PC, ignore, ignore) = processPC(PC_origin, wp_pos);
    MatrixXd c_next = setVertice(c1.transpose(), M_PC).transpose();
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

#endif