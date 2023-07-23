#ifndef HEADER_H
#define HEADER_H

/*
* The global variable for safe polishing 
* Author - Weiye Zhao
* Date - April. 26, 2020
*/

#include <Eigen/Dense>
using namespace Eigen;

// any source file that includes this will be able to use "planner"
extern int planner;
extern MatrixXd pos_idle;
extern MatrixXd pos_middle;
extern MatrixXd theta_ini_polish;
extern MatrixXd theta_ini_measure;
extern MatrixXd Msix2tool;

#endif