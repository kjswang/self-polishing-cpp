#ifndef _MathPack_
#define _MathPack_

/*
* Math operation library for CFS   
* Author - Weiye Zhao
* Date - Nov. 5, 2019
*/

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/LU>
#include "structure.h"
#include "QuadProg++.hh"
// #include "robot_property_eigen.h"
#include <cmath>
#include <vector>
#include <chrono> 
//#include <libqhullcpp/Qhull.h>

// #include "libqhullcpp/RboxPoints.h"
// #include "libqhullcpp/QhullError.h"
// #include "libqhullcpp/QhullQh.h"
// #include "libqhullcpp/QhullFacet.h"
// #include "libqhullcpp/QhullFacetList.h"
// #include "libqhullcpp/QhullFacetSet.h"
// #include "libqhullcpp/QhullLinkedList.h"
// #include "libqhullcpp/QhullPoint.h"
// #include "libqhullcpp/QhullUser.h"
// #include "libqhullcpp/QhullVertex.h"
// #include "libqhullcpp/QhullVertexSet.h"
// #include "libqhullcpp/Qhull.h"

using namespace Eigen;
using namespace std;
//using namespace std::chrono;
// using orgQhull::Qhull;
// using orgQhull::QhullError;
// using orgQhull::QhullFacet;
// using orgQhull::QhullFacetList;
// using orgQhull::QhullFacetListIterator;
// using orgQhull::QhullFacetSet;
// using orgQhull::QhullFacetSetIterator;
// using orgQhull::QhullPoint;
// using orgQhull::QhullPoints;
// using orgQhull::QhullPointsIterator;
// using orgQhull::QhullQh;
// using orgQhull::QhullUser;
// using orgQhull::QhullVertex;
// using orgQhull::QhullVertexList;
// using orgQhull::QhullVertexListIterator;
// using orgQhull::QhullVertexSet;
// using orgQhull::QhullVertexSetIterator;
// using orgQhull::RboxPoints;


void quad2matrix(quadprogpp::Vector<double> u, MatrixXd& u_){
    int row = u.size();
    u_.resize(row,1);
    for(int i=0; i<row; ++i){
        u_(i,0) = u[i];
    }
} 

double norm(MatrixXd& V){
  /*
  * get the norm for a vector represented using Matrix Eigen
  */
    assert(V.cols() == 1 || V.rows() == 1);
    if (V.cols() == 1){
      double sum = 0;
      for (int i=0; i<V.rows(); ++i){
        sum += V(i,0)*V(i,0);
      }
      double sqt = sqrt(sum);
      return sqt;
    }
    if (V.rows() == 1){
      double sum = 0;
      for (int i=0; i<V.cols(); ++i){
        sum += V(0,i)*V(0,i);
      }
      double sqt = sqrt(sum);
      return sqt;
    }
} 


MatrixXd MatPower(MatrixXd& mat, int power){
  MatrixXd new_mat;
  new_mat = mat;
  for (int i=0; i<power-1; ++i){
    new_mat = new_mat * mat;
  }
  if (power == 0)
    new_mat = MatrixXd::Identity(mat.rows(),mat.cols());
  return new_mat;
}

double Matdot(MatrixXd vec1, MatrixXd vec2){
  /*
  * vector dot product of Eigen matrix 
  * mat1 is on the vector 1
  * mat2 is on the vector 2 
  */
  // make sure vec1 and vec2 is actually vector 
  // not some arbitrary sized matrix
  assert(vec1.cols()==vec2.cols() && vec1.rows() == vec2.rows());
  assert(vec1.cols()==1 || vec1.rows()==1);
  // change matrix to vector
  VectorXd v1(Map<VectorXd>(vec1.data(), vec1.cols()*vec1.rows()));
  VectorXd v2(Map<VectorXd>(vec2.data(), vec2.cols()*vec2.rows()));
  double product = v1.dot(v2);
  return product;
}

MatrixXd Matcross(MatrixXd vec1, MatrixXd vec2){
  /*
  * vector cross product of Eigen matrix 
  * mat1 is on the vector 1
  * mat2 is on the vector 2 
  */
  // make sure vec1 and vec2 is actually vector 
  // not some arbitrary sized matrix
  assert(vec1.cols()==vec2.cols() && vec1.rows() == vec2.rows());
  assert(vec1.cols()==1 || vec1.rows()==1);
  // only vector of size 3 is acceptable
  // assert it 
  assert(vec1.cols()*vec2.rows() == 3);
  // change matrix to vector
  Vector3d v1(Map<VectorXd>(vec1.data(), vec1.cols()*vec1.rows()));
  Vector3d v2(Map<VectorXd>(vec2.data(), vec2.cols()*vec2.rows()));
  Vector3d product = v1.cross(v2);
  // transfer eigen vector to eigen matrix
  MatrixXd result(vec1.rows(),vec2.cols());
  result << product(0), product(1), product(2);
  return result;
}



MatrixXd Vcat(MatrixXd mat1, MatrixXd mat2){
  /*
  * vertical concatenation of Eigen matrix 
  * mat1 is on the top 
  * mat2 is on the bottom 
  */
  // if mat1 is empty matrix
  if (mat1.cols() == 0){
    MatrixXd new_mat = mat2;
    return new_mat;
  }
  // mat1 is not empty
  else{
    assert(mat1.cols() == mat2.cols()); // two column number must be same
    MatrixXd new_mat(mat1.rows()+mat2.rows(), mat1.cols()); // <-- D(A.rows() + B.rows(), ...)
    new_mat << mat1, 
               mat2; // <-- syntax is the same for vertical and horizontal concatenation
    return new_mat;
  } 
}


MatrixXd Hcat(MatrixXd mat1, MatrixXd mat2){
  /*
  * horizontal concatenation of Eigen matrix 
  * mat1 is on the left 
  * mat2 is on the right 
  */
  // if mat1 is empty matrix
  if (mat1.rows() == 0){
    MatrixXd new_mat = mat2;
    return new_mat;
  }
  // mat1 is not empty matrix
  else{
    assert(mat1.rows() == mat2.rows()); // two column number must be same
    MatrixXd new_mat(mat1.rows(), mat1.cols()+mat2.cols()); 
    new_mat << mat1, mat2; // <-- syntax is the same for vertical and horizontal concatenation
    return new_mat;
  }
}


double pow_sum(MatrixXd m, int pwr){
  assert(m.cols() == 1);
  double sum = 0;
  MatrixXd tmp;
  tmp = m.array().pow(pwr);
  for (int i=0; i<m.rows(); ++i){
    sum += tmp(i,0);
  }
  return sum;
}


double dot_sum(MatrixXd m, MatrixXd n){
  assert(m.rows() == n.rows() && n.cols()==1 && m.cols()==1);
  double sum = 0;
  for (int i=0; i<m.rows(); ++i){
    sum += m(i,0) * n(i,0);
  }
  return sum;
}


double fixbound(double num){
  if (num < 0){
    num = 0;
  }
  else{
    if (num > 1){
      num = 1;
    }
  }
  return num;
}

void vector2Matrix(MatrixXd &matrix, int row, int col, std::vector<double> v){
  matrix.resize(row, col);
  for (int i=0; i<row; i++){
    for (int j=0; j<col; j++){
      matrix(i,j) = v.at(i*col+j);
    }
  }
}





void setObjValue(quadprogpp::Matrix<double> &G, quadprogpp::Vector<double> &g0, const MatrixXd& HT, const MatrixXd& f){
    int Hn_, fn_;
    Hn_ = HT.rows();
    fn_ = f.rows();
    // cout << "the rows of f is" << fn_ << endl;
    // Hfull
    G.resize(Hn_, Hn_);
    for (int i = 0; i < Hn_; i++) 
        for (int j = 0; j < Hn_; j++)
            G[i][j] = HT(i,j);
    // f

    g0.resize(fn_);
    for (int i = 0; i < fn_; i++){
      // cout << f(i,0) << endl;
      g0[i] = f(i,0);
    }
        
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

    for (int i = 0; i < Sn_; i++){
        ci0[i] = S(i,0);
    }

}

void QPxset(quadprogpp::Vector<double> &x, const MatrixXd& xref){
    x.resize(xref.rows());
    for (int i = 0; i < xref.size(); i++) 
        x[i] = xref(i,0);
}

MatrixXd hatw(MatrixXd w){
  /* hatw is an alternative way of write cross product
  * hatw(w)*a = w crossprodut a 
  */
  assert(w.cols() == 1 && w.rows() == 3);
  double w1 = w(0,0);
  double w2 = w(1,0);
  double w3 = w(2,0);
  MatrixXd hw;
  hw.resize(3,3);
  hw << 0, -1*w3, w2,
        w3, 0, -1*w1,
        -1*w2, w1, 0;
  return hw;
}


bool vecEq(MatrixXd vec1, MatrixXd vec2){
  // judge if two vector is equal 
  MatrixXd diff = vec1 - vec2;
  if (norm(diff) < 0.01)
    return true;
  else
    return false;
}

bool vecNeq(MatrixXd vec1, MatrixXd vec2){
  // judge if two vector is not equal 
  MatrixXd diff = vec1 - vec2;
  if (norm(diff) < 0.01)
    return false;
  else
    return true;
}

MatrixXd cross_point(MatrixXd p1, MatrixXd p2, MatrixXd p3, MatrixXd p4){
    /*
    * get the cross point location of two straight line representd 
    * using two segments p1p2 and p3p4 in 3d space
    */
    // absolute value 
    assert(p1.cols()==1 && p2.cols()==1 && p3.cols()==1 && p4.cols()==1);
    double x1 = p1(0,0);
    double y1 = p1(1,0);
    double z1 = p1(2,0);
    double x2 = p2(0,0);
    double y2 = p2(1,0);
    double z2 = p2(2,0);
    double x3 = p3(0,0);
    double y3 = p3(1,0);
    double z3 = p3(2,0);
    double x4 = p4(0,0);
    double y4 = p4(1,0);
    double z4 = p4(2,0);

    // difference value 
    MatrixXd dp1 = (p1 - p2);
    MatrixXd dp2 = (p3 - p4);
    double dx1 = dp1(0,0);
    double dy1 = dp1(1,0);
    double dz1 = dp1(2,0);
    double dx2 = dp2(0,0);
    double dy2 = dp2(1,0);
    double dz2 = dp2(2,0);

    // analytical solution 
    double x = (dx1*dx2*y3 - dx2*dx1*y1 - dx1*dy2*x3 + dx2*dy1*x1) / (dx2*dy1 - dx1*dy2);
    double y = (dy1*(x - x1))/dx1 + y1;
    double z = (dz1*(y - y1))/dy1 + z1;

    MatrixXd pos(3,1);
    pos << x, y, z;
    return pos;
}

MatrixXd calculateRowMeans(const MatrixXd& matrix) {
    int numRows = matrix.rows();
    int numCols = matrix.cols();

    // Create a column vector to store the row-wise means
    MatrixXd rowMeans(numRows, 1);

    // Calculate the mean for each row
    for (int i = 0; i < numRows; ++i) {
        double sum = 0.0;
        for (int j = 0; j < numCols; ++j) {
            sum += matrix(i, j);
        }
        rowMeans(i, 0) = sum / numCols;
    }

    return rowMeans;
}

MatrixXd calculateColMeans(const MatrixXd& matrix) {
    int numRows = matrix.rows();
    int numCols = matrix.cols();

    // Create a column vector to store the row-wise means
    MatrixXd colMeans(1, numCols);

    // Calculate the mean for each row
    for (int i = 0; i < numCols; ++i) {
        double sum = 0.0;
        for (int j = 0; j < numRows; ++j) {
            sum += matrix(j, i);
        }
        colMeans(0, i) = sum / numRows;
    }

    return colMeans;
}

// Function to calculate the mean of a matrix along a specific dimension
MatrixXd calculateMean(const MatrixXd& mat) {
    return mat.colwise().mean();
}

// Function to apply rotation transformation to a matrix
MatrixXd applyRotation(const MatrixXd& matrix, const Vector3d& axis, double angle) {
    double rad = angle * M_PI / 180.0;
    AngleAxisd rotation(rad, axis);
    MatrixXd rotatedMatrix = rotation.toRotationMatrix() * matrix;
    return rotatedMatrix;
}

// Function to apply translation transformation to a matrix
MatrixXd applyTranslation(const MatrixXd& matrix, const Vector3d& translation) {
    Matrix4d translationMatrix = Matrix4d::Identity();
    translationMatrix.block<3, 1>(0, 3) = translation;
    MatrixXd translatedMatrix = translationMatrix * matrix;
    return translatedMatrix;
}

MatrixXd mtx_translate(MatrixXd t){
  MatrixXd T_t = MatrixXd::Identity(4,4);
  if (t.rows() == 3 && t.cols() == 1) {
        for (int i = 0; i < 3; ++i) {
            T_t(i, 3) = t(i, 0);
        }
    } else {
        std::cerr << "Error: 't' should be a 3x1 matrix (column vector)." << std::endl;
    }
    return T_t;
}

MatrixXd mtx_rotate(MatrixXd n, MatrixXd a){
    double theta = n.row(0).norm();
    MatrixXd T_t_a = mtx_translate(a);

    MatrixXd k = n;
    if (theta != 0){
      k = n/theta;
    }

    double k_x = k(0,0);
    double k_y = k(0,1);
    double k_z = k(0,2);

    MatrixXd K(3,3);
    K << 0, -k_z, k_y,
        k_z, 0, -k_x,
        -k_y, k_x, 0;

    MatrixXd R = MatrixXd::Identity(3,3) + sin(theta)*K + (1-cos(theta))*K*K;
    MatrixXd T_r = MatrixXd::Identity(4,4);
    T_r.block(0,0,3,3) = R;
    T_r = mtx_translate(a) * T_r * mtx_translate(-a);
    return T_r;
}


#endif