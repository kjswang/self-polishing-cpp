#ifndef _DISCON3D_H_
#define _DISCON3D_H_

/*
* Robot-obstacle distance calculation library 
* Author - 
* Date - 
*/
#include<iostream>
#include <vector>
#include <Eigen/Dense>
#include "structure.h"
#include "CapPos.h"
#include <cmath>
#include "math_pack.h"
#include <libqhullcpp/Qhull.h>
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"
using namespace std;
using namespace Eigen;
using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetSetIterator;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexSetIterator;
using orgQhull::RboxPoints;

double distLinSeg(MatrixXd s1, MatrixXd e1, MatrixXd s2, MatrixXd e2){
  MatrixXd d1, d2, d12;
  d1 = e1 - s1;
  d2 = e2 - s2;
  d12 = s2 - s1;
  
  double D1 = pow_sum(d1, 2);
  double D2 = pow_sum(d2, 2);
  double S1 = dot_sum(d1, d12);
  double S2 = dot_sum(d2, d12);
  double R = dot_sum(d1, d2);
  double den = D1*D2 - pow(R,2);


  //if one of the segments is a point
  double u,t,uf;
  if (D1 == 0 || D2 == 0){
    if (D1 != 0){
      u = 0;
      t = S1 / D1;
      t = fixbound(t);
    }
    else if (D2 != 0){ 
      t = 0;
      u = -S2 / D2;
      u = fixbound(u);
    }
    else{ // both line segment is point 
      t = 0;
      u = 0;
    }
  }

  else if (den == 0){ // line is parallel 
    t = 0;
    u = -S2/D2;
    uf = fixbound(u);
    if (uf != u){
      t = (uf*R+S1)/D1;
      t = fixbound(t);
      u = uf;
    }
  }

  else{ // general case
    t = (S1*D2-S2*R)/den;
    t = fixbound(t);
    u = (t*R-S2)/D2;
    uf = fixbound(u);
    if (uf != u){
      t = (uf*R+S1)/D1;
      t = fixbound(t);
      u = uf;
    }
  }

    // compute the distance, given t and u 
    MatrixXd distm;
    distm = d1*t-d2*u-d12;
    double dist = distm.norm();
    // point
    MatrixXd point(3,2);
    point.block(0,0,3,1) = s1 + d1*t;
    point.block(0,1,3,1) = s2 + d2*u;
    return dist;

}


double dist_arm_3D_Heu(MatrixXd& theta, MatrixXd& DH, MatrixXd& base, MatrixXd& obs, capsule cap[]){
  int nstate = theta.rows();
  for (int i=0; i<nstate; ++i){
    DH(i,0) = theta(i,0);
  }
  double d = std::numeric_limits<double>::infinity(); //set d as positive infinity

  // forward kinematics for position calculation under theta
  lineseg pos[DH.rows()];
  MatrixXd M[DH.rows()+1];
  CapPos(base, DH, cap, M, pos);
  // cout << base << endl << DH << endl;
  
  // calculate the closest distance 
  double dis;
  int link_id;
  for (int i=0; i<nstate; ++i){
    dis = distLinSeg(pos[i].p1, pos[i].p2, obs.block(0,0,3,1), obs.block(0,1,3,1)); // this is 3d scenario. where obstacle position in 3d form
    if (dis < d){
      d = dis;
      link_id = i+1;
    }
  }
  return d;
}


double distLinSeg2PC(MatrixXd P, MatrixXd Q, MatrixXd PC){
    MatrixXd PQ = P - Q;
    MatrixXd PC_Q = PC.colwise() - Q.col(0);
    if (PQ.norm() == 0) {
        return PC_Q.colwise().norm().minCoeff();
    }

    MatrixXd t = (PC_Q.transpose()*PQ) / (PQ.transpose()*PQ).norm();
    t = t.cwiseMax(0);
    t = t.cwiseMin(1);

    MatrixXd TQ = PQ*t.transpose();
    MatrixXd dist = (PC_Q - TQ).colwise().norm();

    // MatrixXd T = PC_Q.transpose() * PQ / (PQ.transpose() * PQ);
    // T = T.cwiseMax(0).cwiseMin(1);
    // MatrixXd TQ = PQ.transpose() * T.transpose();
    // MatrixXd dist = (PC_Q - TQ).colwise().norm();
    return dist.minCoeff();
}

double distLinSeg2PC_arma(mat P, mat Q, mat PC) {
    arma::mat PQ = P - Q;
    arma::mat PC_Q = PC.each_col() - Q;

    if (arma::norm(PQ) == 0) {
        arma::mat norms = arma::sqrt(arma::sum(arma::square(PC_Q), 0));
      return norms.min();
    }

    arma::mat T = (arma::trans(PC_Q) * PQ) / arma::norm(arma::trans(PQ) * PQ);
    T.elem(arma::find(T < 0)).zeros();
    T.elem(arma::find(T > 1)).ones();

    arma::mat norms = arma::sqrt(arma::sum(arma::square(PC_Q - PQ*T.t()), 0));
    return norms.min();
}


double dist_arm_PC(MatrixXd theta, MatrixXd DH, MatrixXd base, capsule cap[], tool tl, MatrixXd PC, int check = 0){
    
    int nlink = DH.rows();
    DH.block(0,0,nlink,1) = theta;
    double d = INFINITY;

    // lineseg* pos = CapPos_real(base, DH, cap, tl);
    PC.transposeInPlace();

    lineseg* pos = CapPos_real(base, DH, cap, tl);

    //CapPos(base, DH, cap, M, pos);
    MatrixXd point1, point2;
    double dist, radius;
    vector<double> distList(sizeof(pos));
    for (int i = 0; i < sizeof(pos); i++) {
        point1 = pos[i].p1;
        point2 = pos[i].p2;
        radius = pos[i].r;
        dist = distLinSeg2PC(point1, point2, PC);
        dist -= radius;
        distList[i] = dist;
    }
    d = *std::min_element(distList.begin(),distList.end());
    return d;
}


MatrixXd convhulln(MatrixXd points){
    vector<double> vectorData;
    for (int i = 0; i < points.rows(); ++i) {
        for (int j = 0; j < points.cols(); ++j) {
            vectorData.push_back(points(i, j));
        }
    }


    Qhull myqhull;
    // myqhull.runQhull("m", 3, points.size() / 3, points.data(), "Qt Qa Qc Qi Qw Qs Q0 Q3 Q4 Q6 Q12");
    myqhull.runQhull("m", 3, vectorData.size() / 3, vectorData.data(), "Qt");
    QhullFacetList facets= myqhull.facetList();
    std::vector<std::vector<int> > facetVertices;
    cout << "num_facets" << myqhull.facetCount() << endl;
    // for(QhullFacet f : facets)
    QhullFacetListIterator j(facets);
    while(j.hasNext()){
        QhullFacet f= j.next();
        std::vector<int> vertices;
        if(!f.isGood()){
            // ignore facet
        }else if(!f.isTopOrient() && f.isSimplicial()){ /* orient the vertices like option 'o' */
            QhullVertexSet vs= f.vertices();
            vertices.push_back(vs[1].point().id());
            vertices.push_back(vs[0].point().id());
            for(int i= 2; i<(int)vs.size(); ++i){
                vertices.push_back(vs[i].point().id());
            }
            facetVertices.push_back(vertices);
        }else{  /* note: for non-simplicial facets, this code does not duplicate option 'o', see qh_facet3vertex and qh_printfacetNvertex_nonsimplicial */
            // for(QhullVertex vertex : f.vertices()){

            QhullVertexSetIterator k(f.vertices());
            while(k.hasNext()){
                QhullVertex vertex= k.next();
                QhullPoint p= vertex.point();
                vertices.push_back(p.id());
            }
            facetVertices.push_back(vertices);
        }
    }
    int numRows = facetVertices.size();
    Eigen::MatrixXd matrix(numRows, 3);

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < 3; ++j) {
            matrix(i, j) = facetVertices[i][j];
        }
    }
    return matrix;
}

bool inhull(const Eigen::MatrixXd &testpoints, const Eigen::MatrixXd &hullpoints) {

    MatrixXd tess = convhulln(hullpoints);
    cout << tess << endl;
    int nt, c;
    nt = tess.rows();
    c = tess.cols();
    double tol = 0.0;
    MatrixXd ab(nt, 3);

    // Eigen::MatrixXd tessIndices1 = hullpoints.row(tess.col(0)).array();  // Extract rows using indices
    // Eigen::MatrixXd tessIndices2 = hullpoints.row(tess.col(1)).array();

    // ab = tessIndices1 - tessIndices2;
    for (int i=0; i<nt; i++){
      int index1 = tess(i,0);
      int index2 = tess(i,1);

      ab.row(i) = hullpoints.row(index1) - hullpoints.row(index2);
    }

    MatrixXd ac(nt,3);

    for (int i = 0; i < nt; ++i) {
        int index3 = tess(i, 0);
        int index4 = tess(i, 2);

        ac.row(i) = hullpoints.row(index3) - hullpoints.row(index4);
    }

    MatrixXd nrmls(nt,3);
    for (int i=0; i<nt;i++){
      Vector3d v = ab.row(i);
      Vector3d w = ac.row(i);
      nrmls.row(i) = v.cross(w);
    }

    MatrixXd degenflag(nt,1);
    degenflag.setZero();
    MatrixXd nrmllen = nrmls.rowwise().norm().col(0);
    MatrixXd nrmllen_reci = nrmllen.array().inverse();
    MatrixXd nrmllen_new = nrmllen_reci.replicate(1,3);
    MatrixXd vec_norm = nrmls.array() * nrmllen_new.array();

    MatrixXd center = hullpoints.colwise().mean().row(0);

    MatrixXd a(nt,3);
    for (int i = 0; i < nt; ++i) {
        int index5 = tess(i, 0);
        a.row(i) = hullpoints.row(index5);
    }

    MatrixXd rep_center = center.replicate(nt,1);

    MatrixXd dp = ((rep_center-a).array() * vec_norm.array()).rowwise().sum();

    MatrixXi k = (dp.array() < 0).cast<int>();

    for (int i=0; i < k.rows(); i++){
      if (k(i,0) == 1){
        vec_norm.row(i) = -vec_norm.row(i);
      }
    }

    MatrixXd aN = (vec_norm.array() * a.array()).rowwise().sum();

    int in = 0;
    double memblock = 1e6;
    int blocks = max(1, static_cast<int>(floor(1/(memblock/nt))));
    MatrixXd aNr = aN;
    MatrixXd temp = vec_norm*testpoints.transpose() - aNr;

    Eigen::Matrix<bool, Eigen::Dynamic, 1> condition_met = (temp.array() >= -tol).colwise().all();
    return condition_met(0,0);
}


#endif