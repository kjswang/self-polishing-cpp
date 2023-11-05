#ifndef HEADER_H
#define HEADER_H

/*
* The global variable for safe polishing 
* Author - Weiye Zhao
* Date - April. 26, 2020
*/

#include <Eigen/Dense>
#include <armadillo>
using namespace Eigen;
using namespace arma;

// any source file that includes this will be able to use "planner"
extern int planner;
extern MatrixXd pos_idle;
extern MatrixXd pos_middle;
extern MatrixXd theta_ini_polish;
extern MatrixXd theta_ini_measure;
extern MatrixXd Msix2tool;

extern mat theta_ini_polish_arma;
extern mat theta_ini_measure_arma;
extern mat Msix2tool_arma;

extern mat fmincon_x;
extern int ft;

#endif