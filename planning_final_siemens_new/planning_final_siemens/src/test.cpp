/*
 * developer testing code
 * just ignore it.
 */

#include <iostream>
#include <vector>
// #include <Eigen/Dense>
#include <Eigen/Geometry> 
#include <chrono> 
#include <string>
#include <map>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

int main() {

  // MatrixXd A, B, C;
  // A = MatrixXd::Identity(1000,1000);
  // B = MatrixXd::Identity(1000,1000);
  // auto start_obj = high_resolution_clock::now();
  // C = A * B;
  // auto end_obj = high_resolution_clock::now();
  // auto duration_obj = duration_cast<microseconds>(end_obj - start_obj); 
  // cout << "the procedure time is: " << double(duration_obj.count())/1000000.0 << "seconds" << endl; 
  // return 0;

  Matrix3d mat;
  mat = MatrixXd::Identity(3,3);
  Quaterniond q;
  q = mat;
  cout << "the x is" << q.x() << endl;
  cout << "the y is" << q.y() << endl;
  cout << "the z is" << q.z() << endl;
  cout << "the w is" << q.w() << endl;
  Matrix3d mat1 = mat*-1;
  q = mat1;
  cout << "the x is" << q.x() << endl;
  cout << "the y is" << q.y() << endl;
  cout << "the z is" << q.z() << endl;
  cout << "the w is" << q.w() << endl;
  q = mat;
  cout << "the x is" << q.x() << endl;
  cout << "the y is" << q.y() << endl;
  cout << "the z is" << q.z() << endl;
  cout << "the w is" << q.w() << endl;
  // MatrixXd rot;
  // rot = MatrixXd::Identity(3,3);
  // // Quaterniond q(rot);
  // Quaterniond q = Quaterniond(rot);
  // // q = rot;
  // cout << q.x();

  // map<int, int> gquiz1; 
  // gquiz1.insert(pair<int, int>(1, 40)); 
  // gquiz1.insert(pair<int, int>(2, 30)); 
  // gquiz1.insert(pair<int, int>(3, 60)); 

  // map<int, int>::iterator itr; 
  // for (itr = gquiz1.begin(); itr != gquiz1.end(); ++itr) { 
  //       cout << '\t' << itr->first 
  //            << '\t' << itr->second << '\n'; 
  //   } 
  //   cout << endl;
 

}



  