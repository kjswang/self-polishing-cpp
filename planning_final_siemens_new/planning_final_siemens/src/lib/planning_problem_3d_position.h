#ifndef _PLANNINGPROBLEM3D_H_
#define _PLANNINGPROBLEM3D_H_

/*
* Planning problem class defintion 
* CFS solving joint position reference
* Author - Weiye Zhao
* Date - Nov. 5, 2019
*/

#include "customized_print.h"
#include "load_parameter.h"
#include <Eigen/Dense>
#include "robot_property_eigen.h"
#include "distance_constraint_3d.h"
#include "structure.h"
#include "numerical_gradient.h"
#include <vector>
#include "math_pack.h"
using namespace std;

// Motion Planning Problem
class PlanningProblem {
public:

    int nstep_;
    int dim_; // dim = nstate
    int nobs_;
    int neq_;
    int njoint_;
    int nstate_;
    int nu_;
    int horizon_;
    double margin_;
    // MatrixXd DH_; 
    MatrixXd obs_;
    MatrixXd xref_;
    MatrixXd xori_;
    MatrixXd base_;

    //Kinematic matrix
    MatrixXd Baug_;
    vector<double> weight_ref_;
    vector<double> weight_self_;
    

    PlanningProblem(Robot& robot){
        // setting parameter, and reference.
        nobs_ = 1; // current set obstacle number to 1
        vector<double> xref_vec;
        loadXrefPosition(xref_vec);
        setReferenceTrajectory(xref_vec, xref_);
        base_.resize(3,1);
        loadPlanningProblemSetting(njoint_, nstate_, nu_, margin_, neq_, base_, weight_ref_, weight_self_);
        dim_ = njoint_; // dimension
        obs_.resize(3,2);
        obs_ << 4.106, 4.106,
                8.313, 8.313,
                1.072, 1.538;
        horizon_ = xref_.rows() / dim_; // planning horizon
        nstep_ = horizon_;
        
        // setting xori_
        xori_ = xref_;
        
    }


    void setReferenceTrajectory(vector<double>& ref_vec, MatrixXd& ref){
        ref.resize(ref_vec.size(),1);
        for (int i=0; i<ref_vec.size(); ++i){
            ref(i,0) = ref_vec[i];
        }
    }


    int setCostMatrix3d(vector< vector<double> >& Qref, vector< vector<double> >& Qself){
        vector< vector<double> > Pos_(Qref);
        vector< vector<double> > Vel_(Qref);
        vector< vector<double> > Adiff_((nstep_-2)*dim_, vector<double>(nstep_*dim_));
        vector< vector<double> > Acc_(Qref);

        for (int i=0; i < nstep_*dim_; i++)
        {
            Pos_[i][i] = 1;
            if ((i >= dim_) && (i<(nstep_-1)*dim_)){
            Vel_[i][i] = 2;
            Vel_[i][i-dim_] = -1;
            Vel_[i-dim_][i] = -1;
            Vel_[i][i+dim_] = -1;
            Vel_[i+dim_][i] = -1;
            }else{
                Vel_[i][i] = 1;
            }
        }
        cout << "pos and vel matrix set" << endl;

        // auto start_obj = high_resolution_clock::now();

        for (int i=0; i < nstep_*dim_; i++)
        { 
            if (i < dim_){
                Acc_[i][i] = 1;
                Acc_[i][i+dim_] = -2;
                Acc_[i][i+dim_*2] = 1;
                Acc_[i+dim_][i] = -2;
                Acc_[i+dim_*2][i] = 1;
                continue;
            }

            if(i >= dim_ && i < dim_*2){
                Acc_[i][i] = 5;
                Acc_[i][i+dim_] = -4;
                Acc_[i][i+dim_*2] = 1;
                Acc_[i+dim_][i] = -4;
                Acc_[i+dim_*2][i] = 1;
                continue;
            }

            if(i >= dim_*2 && i < dim_*(nstep_-2)){
                Acc_[i][i] = 6;
                Acc_[i][i+dim_] = -4;
                Acc_[i][i+dim_*2] = 1;
                Acc_[i+dim_][i] = -4;
                Acc_[i+dim_*2][i] = 1;
                continue;
            }

            if(i >= dim_*(nstep_-2) && i < dim_*(nstep_-1)){
                Acc_[i][i] = 5;
                Acc_[i][i+dim_] = -2;
                Acc_[i+dim_][i] = -2;
                continue;
            }

            if(i >= dim_*(nstep_-1)){
                Acc_[i][i] = 1;
            }
        }
        cout << "Acc matrix set" << endl;

        for (int i=0; i < nstep_*dim_; i++){
            for (int j=0; j < nstep_*dim_; j++){
                Qref[i][j] = weight_ref_[0] * Pos_[i][j] / nstep_ + weight_ref_[1] * Vel_[i][j] * (nstep_-1)  + weight_ref_[2] * Acc_[i][j] * pow(nstep_-1,4) / (nstep_ - 2);
                Qself[i][j] = weight_self_[0] * Pos_[i][j] / nstep_ + weight_self_[1] * Vel_[i][j] * (nstep_-1) + weight_self_[2] * Acc_[i][j] * pow(nstep_-1,4) / (nstep_ - 2);
            }
        }
        return 0;
    }


    void setObjective3d(MatrixXd& H, MatrixXd& f, Robot& robot){
        vector< vector<double> > H_(nstep_*dim_, vector<double>(nstep_*dim_));
        vector<double> f_(nstep_*(dim_+neq_));
        vector< vector<double> > Qref(nstep_*dim_,vector<double>(nstep_*dim_));
        vector< vector<double> > Qself(nstep_*dim_,vector<double>(nstep_*dim_));

        setCostMatrix3d(Qref, Qself); 
        for (int i=0; i<nstep_*dim_; i++){
            f_[i] = 0;
            for (int j=0; j<nstep_*dim_; j++)
            {
                H_[i][j] = Qref[i][j] + Qself[i][j];
                f_[i] += -1*Qref[i][j]*xref_(j,0);
            }
        }
     
        //transform Vector-matrix into matrix
        H.resize(nstep_ * dim_, nstep_ * dim_);
        for (int i=0; i<nstep_*dim_; ++i){
            for (int j=0; j<nstep_*dim_; ++j){
                H(i,j) = H_[i][j];
            }
        }
        //transform Vector-vector into matrix
        // f.resize(nstep_*dim_,1);
        f.resize(1, nstep_*dim_);
        for (int i=0; i<nstep_*dim_; ++i){
            f(0,i) = f_[i];
        }
    }


    /*
    * set the linear inequality constriant for CFS QP solver 
    */
    int linConstraint(const Robot& robot, MatrixXd& Lfull, MatrixXd& S){
        Lfull.resize(horizon_, horizon_*dim_);
        S.resize(horizon_, 1);
        S = MatrixXd::Zero(horizon_, 1);
        // create new variable space 
        MatrixXd theta, DH_use;
        capsule cap[6];
        for (int i=0; i<6; ++i){
            cap[i] = robot.cap[i];
        }
        
        // start constraint specification
        MatrixXd grad, s, l, d, Diff;
        double distance;
        for (int i=0; i<horizon_; ++i){
        // for (int i=0; i<1; ++i){
            theta = xref_.block(i*dim_, 0, dim_, 1);
            DH_use = robot.DH.block(0,0,dim_,robot.DH.cols());
            distance = dist_arm_3D_Heu(theta, DH_use, base_, obs_, cap);
            d.resize(1,1);
            d(0,0) = distance - margin_;
            // numerical gradient
            grad.resize(1, theta.rows());
            grad = central_diff(theta, DH_use, base_, obs_, cap);
            Diff.resize(1, horizon_*dim_); // construct gradient of xref in terms of every time step
            Diff = MatrixXd::Zero(1, horizon_*dim_); // the rest is zero
            Diff.block(0, i*dim_, 1, dim_) = grad;

            // reset linear constriants 
            s = d - Diff*xref_;
            l = -1*Diff;

            // concatnate s, l to construct S, Lfull
            assert(l.cols() == horizon_*dim_);
            assert(s.cols() == 1 && s.rows() == 1);
            S(i,0) = s(0,0);
            Lfull.block(i, 0, 1, l.cols()) = l;
        }
      return 0;
    }


    /*
    * set the linear equality constriant for CFS QP solver 
    */
    void eqConstraint(MatrixXd& Aeq, MatrixXd& beq, int fix=1){
        int count = 0;
        Aeq.resize(dim_*2, nstep_*dim_);
        beq.resize(dim_*2,1);
        Aeq = MatrixXd::Zero(dim_*2, nstep_*dim_);
        beq = MatrixXd::Zero(dim_*2, 1);
        for (int i=0; i<dim_ * fix; i++){
            Aeq(count,i) = 1;
            beq(count,0) = xref_(i,0);
            count++;
            Aeq(count, (nstep_-fix)*dim_+i) = 1;
            beq(count, 0) = xref_((nstep_-fix)*dim_+i,0);
            count++;
        }
    }


    double getCost(const MatrixXd& H_, const MatrixXd& f_){
        double cost = 0;
        MatrixXd obj;
        obj = 0.5 * xref_.transpose() * H_ * xref_ + f_ * xref_;
        assert(obj.cols()==1 && obj.rows()==1);
        cost = obj(0,0);
        return cost;
    }
};
#endif