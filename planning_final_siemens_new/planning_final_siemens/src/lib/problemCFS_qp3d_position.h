#ifndef _PROBLEMCFS3DPOS_H_
#define _PROBLEMCFS3DPOS_H_
/*
* CFS class defintion 
* CFS solving joint position reference
* Author - Weiye Zhao
* Date - Nov. 5, 2019
*/
#pragma once
#include "customized_print.h"
#include "PlanningProblem_3d_position.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include "robot_property_eigen.h"
#include "math_pack.h"

using namespace std;

class ProblemCFS {


    public:
        PlanningProblem* pp_;  
        double iteration_time_;
        int iteration_;
        double tmp_cost;
        vector<double> qp_time_;
        vector<double> process_time_;
        vector<double> cost_;
        MatrixXd soln_;
        MatrixXd xold_;
        MatrixXd xR_; // xref is the same as xR, by stacking each column of xR;

    ProblemCFS(PlanningProblem* pp) {
        pp_ = pp;
    }

    int iteration(int iterMax, double tolX, Robot& robot){

        // set Objective from planning problem
        MatrixXd Hfull_, f_;
        MatrixXd Lfull_, S_; // inequality constraints
        MatrixXd Aeq_, beq_; // equality constraints
        MatrixXd nLT_, nAeqT_; // negative transpose matrix
        quadprogpp::Matrix<double> G, CE, CI;
        quadprogpp::Vector<double> g0, ce0, ci0, x;

        // initialize reference input
        QPxset(x, pp_->xref_);// set initial value of u;
        xold_ = pp_->xref_;

        // set objecive of QP solver 
        // auto start_obj = high_resolution_clock::now();
        pp_->setObjective3d(Hfull_,f_,robot);
        // auto stop_obj = high_resolution_clock::now(); 
        // auto duration_obj = duration_cast<microseconds>(stop_obj - start_obj); 
        // cout << "the obj time is: " << double(duration_obj.count())/1000000.0 << "seconds" << endl; 
        setObjValue(G, g0, Hfull_, f_);

        // record time
        double iter_start_time = 0;
        double iter_end_time = 0;
        double qp_start_time = 0;
        double qp_end_time = 0;
        double process_start_time = 0;
        double process_end_time = 0;

        // The iteration start
        qp_time_.clear();
        cost_.clear();
        iter_start_time = clock();     // record time

        for (int k=0; k<iterMax; k++){
            // for (int k=0; k<1; k++){
            cout << "----------------------" << endl;
            cout << "Iteration " << k << endl;

            // Processing
            // reset the QP Quadratic term and linear term
            // mandatory: G should be reset after solving one QP problem
            setObjValue(G, g0, Hfull_, f_);

            // set linear inequality constraint
            process_start_time = clock(); 
            pp_->linConstraint(robot, Lfull_, S_);
            nLT_ = -1*Lfull_.transpose(); // negative transpose
            setConstraint(CI, ci0, nLT_, S_); 

            // set linear equality constraint 
            pp_->eqConstraint(Aeq_, beq_);
            nAeqT_ = -1*Aeq_.transpose();
            setConstraint(CE, ce0, nAeqT_, beq_); 
            /* uncomment blow to cancel equality constraints*/
            // CE.resize(pp_->horizon_*pp_->njoint_,0);
            // ce0.resize(0);
            // clock stop
            process_end_time = clock();
            process_time_.push_back((process_end_time - process_start_time)/CLOCKS_PER_SEC);


            // Solve the subproblem
            printMessage("solving QP...");
            qp_start_time = clock();
            tmp_cost = solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
            quad2matrix(x, pp_->xref_); // update uref_ in planning problem;
            // clock stop
            qp_end_time = clock();
            qp_time_.push_back((qp_end_time - qp_start_time)/CLOCKS_PER_SEC);
            cout << "cost temporal cost is: " << tmp_cost << endl;
            cost_.push_back(tmp_cost);


            //Convergence check
            // compare new xreference
            if (checkConvergence(pp_->xref_,xold_,tolX)){
                cout << "Converged at step " << k+1 << endl;
                iteration_ = k;
                break;
            }
            // check if inequality constraints are satisifed
            // if (checkLocalOpt(Lfull_, S_, pp_->xref_)){
            //   cout << "Local Optima found at step " << k << endl;
            //   iteration_ = k;
            //   break;
            // }
            // record the current reference trajectory
            xold_ = pp_->xref_;
        }
        soln_ = pp_->xref_;
        iter_end_time = clock();
        iteration_time_ = (iter_end_time - iter_start_time)/CLOCKS_PER_SEC;
        return 0;
    }


    bool checkConvergence(const MatrixXd& x, MatrixXd& xold, const double& tolX){
        MatrixXd diff;
        diff = x - xold;
        double norm_diff;
        norm_diff = diff.norm();
        if (norm_diff < tolX){
            return true;
        }
        else{
            return false;
        }
    }


    bool checkLocalOpt(const MatrixXd& L, const MatrixXd& S, const MatrixXd& x){
        bool optima = true;
        MatrixXd diff;
        diff = S - L*x;
        for (int i=0; i<diff.rows(); ++i){
            if (diff(i,0) <= 0){
                optima = false;
                break;
            }
        }
        return optima;
    }


    void QPxset(quadprogpp::Vector<double> &x, const MatrixXd& xref){
        x.resize(xref.rows());
        for (int i = 0; i < xref.size(); i++) 
            x[i] = xref(i,0);
    }


    void setConstraint(quadprogpp::Matrix<double> &CI, quadprogpp::Vector<double> &ci0, const MatrixXd& LT, const MatrixXd& S){
    /* push Lfull_ to qp solver
       be careful that:
       according to QP solver, CI should be -Lfull_*/
        int Lnr_ = LT.rows();
        int Lnc_ = LT.cols();
        int Sn_ = S.rows();
        CI.resize(Lnr_, Lnc_);
        for (int i = 0; i < Lnr_; i++) 
            for (int j = 0; j < Lnc_; j++)
                CI[i][j] = LT(i,j);

        // f
        ci0.resize(Sn_);
        for (int i = 0; i < Sn_; i++) 
            ci0[i] = S(i,0);
    }


  void printResult(){
    cout << "the process time (second) is:";
    printVector(process_time_.data(),process_time_.size());
    cout << "the QP time (second) is:";
    printVector(qp_time_.data(),qp_time_.size());
    // cout << "the costs are:";
    // printVector(cost_.data(),cost_.size());
    cout << "total iteration time = " << iteration_time_ << "second" << endl;
  }
};

#endif
