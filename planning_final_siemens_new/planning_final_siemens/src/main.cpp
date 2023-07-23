/*
* The ROS node publish the safe polish reference trajectory
* include libraries of safe polishing library
* Author - Jiansong Wang
* Date - 
*/

// #include "ros/ros.h"
// #include <nav_msgs/Path.h>

#include <iostream>
#include <Eigen/Dense>

#include "safe_polish_HTPC.h"
#include "CapPos.h"
#include "load_param.h"
#include <fstream>
#include "structure.h"
#include "global_var.h"
const static IOFormat txtFormat(StreamPrecision, DontAlignCols, ",", "\n");

using namespace std;
using namespace Eigen;

MatrixXd theta_init_default(6,1);
MatrixXd wp_pos_init_default(3,1);
MatrixXd start2exe_traj, execution_traj, exit_traj;
int main(){

    cout << "beforewestart" << endl;

    safe_polish(start2exe_traj, execution_traj, exit_traj);


    return 0;
}
