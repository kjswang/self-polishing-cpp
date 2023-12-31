/*
* The ROS node publish the safe polish reference trajectory
* include libraries of safe polishing library
* Author - Jiansong Wang
* Date - 
*/


//new

// #include "ros/ros.h"
// #include <nav_msgs/Path.h>

#include <iostream>
#include <Eigen/Dense>

#include "safe_polishing.h"
#include "CapPos.h"
#include "load_parameter.h"
#include <fstream>
#include "structure.h"
#include "global_var.h"
#include <libqhullcpp/Qhull.h>
#include <vector>


const static IOFormat txtFormat(StreamPrecision, DontAlignCols, ",", "\n");

using namespace std;
using namespace Eigen;

MatrixXd theta_init_default(6,1);
MatrixXd wp_pos_init_default(3,1);
MatrixXd start2exe_traj, execution_traj, exit_traj;

int planner; // 0 is for measurement and 1 is for polishing
MatrixXd pos_idle(6,1); // the idle pos
MatrixXd pos_middle(6,1); // the intermediate pos 
MatrixXd theta_ini_polish(6,1);
MatrixXd theta_ini_measure(6,1);
MatrixXd Msix2tool(4,4);


int main(int argc, char **argv){

    loadjnt2tool(Msix2tool);
    planner = 1;
    
    safe_polish(start2exe_traj, execution_traj, exit_traj);


    return 0;
}
