/*
* CFS main function 
* include libraries of planning problem class, problem CFS class, etc. 
* CFS demonstration, performance recorded, results saved. 
* Author - Weiye Zhao
* Date - Nov. 10, 2019
*/

#include "PlanningProblem_3d_position.h"
#include "ProblemCFS_qp3d_position.h"
#include "loadParameter.h"
#include "customizedPrint.h"
#include <cstdio>
#include <ctime>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include "QuadProg++.hh"
#include <fstream>
#include <iterator>
#include <vector>
#include <array>
#include "robot_property_eigen.h"
#include <Eigen/Dense>
#include <string>

using namespace std;
using namespace Eigen;
/**
 * CFS main function.
 */

int main() {
    /*
    * robot property initialization, use robot type number specification
    * 1 - M16iB
    * 2 - GP50
    */
    string name = "M16iB"; 
    Robot robot(name);

    // CFS planning problem initialization
    PlanningProblem pp(robot);
    ProblemCFS cfs(&pp);
    cout << "solving" << endl;
    // start execution 
    double tic = clock();
    cfs.iteration(loadMaxIt(), 0.01, robot);
    double toc = clock();

    // print results
    cout << "the plannig horizon is: " << pp.nstep_ << endl;
    // cout << "the size of solution is: " << cfs.soln_.rows() << endl;
    cfs.printResult();

    /*
    Save the reference trajectory solution
    */
    std::ofstream ofs;
    ofs.open("xreference.txt", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    std::ofstream outfile;
    outfile.open("xreference.txt", std::ios_base::app);
    outfile << "\n";
    for (int i=0; i<cfs.soln_.rows(); ++i){
        outfile << cfs.soln_(i,0);
        outfile << "\n";
    }
}


