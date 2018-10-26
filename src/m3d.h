/* m3d
 * 	These files may be found at:
 * 	http://software.seg.org/2017/0006
 *
 * 	Please first read and accept:
 * 	http://software.seg.org/disclaimer.txt
*/

/*
  Chen Huang,    csuchenhuang@csu.edu.cn
  Zhengyong Ren, zhengyong@csu.edu.cn
  Central South University of China, 2018
*/

#ifndef  _M3D_H
#define  _M3D_H

// C++ includes
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "node.h"
#include "tet.h"
#include "grid.h"

/*Revised introduction: 
1) rewrited the compute_1_R3 function.
2) rewrited the compute_grad_grad_1/R function.
3) added compute_xpxq_R5\compute_xq_R5\compute_1_R5 surface integral function
   and line_1_R3\line_mij_r2R\line_mij_r2R3 three line integral function.*/

class M3D 
{
 public:
  M3D(Tet& tet, Point r);
  void V(double& v);
  void B(Point& b);
  void T(double tb[3][3]);
  double sgn(const double m);
  void set(double a[3][3], double b);
  void compute_1_R(unsigned int i, double &k0); 
  void compute_grad_1_R(unsigned int i, Point &k1); 
  void compute_grad_grad_1_R(unsigned int i, double k2[3][3]);
  void compute_1_R3(unsigned int i,  double &k3);
  void compute_xpxq_R5(unsigned int i, unsigned int p, unsigned int q, double &k4);
  void compute_xq_R5(unsigned int i, unsigned int q, double &k5);
  void compute_1_R5(unsigned int i, double &k6);
  void line_1_R3(unsigned int i, unsigned int j, double &k7);
  void line_mij_r2R(unsigned int i, unsigned int j, double &k8);
  void line_mij_r2R3(unsigned int i, unsigned int j, double &k9);
  double get_tri_size(std::vector<Point>& T);
  const double TOL;
  const double Cm ; 

 private:
  // variables for edges 
  double    Rij0     [4][3]; 
  double    Rij1     [4][3];  
  double    Sij0     [4][3];  
  double    Sij1     [4][3];  
  double    Rij      [4][3];  
  double    mij      [4][3]; // mij is constant along edge eij
  Point     mij_hat  [4][3]; // mij_hat is the outgoing normal vector of edge eij
  Point     eij_hat  [4][3]; // eij_hat is the tangential normal vector of edge eij
  // adding two edge variables
  Point     rij_n    [4][3]; // rij_n is the projection point of the site r on edge Cij
  Point     rhoij_n  [4][3]; // is the unit vector from oi to point rij_n   
  // variables for facets
  Point     ni_hat   [4];    // ni_hat is the outgoing normal vector of surface Ti
  double    hi       [4];    // equation (8)
  Point     oi       [4];    // oi is the projection point of site r on surface Ti
  // observation site and source element
  Tet&               tet;    // a tetrahedron
  Point              r;      // an observation site, which is denoted by r' in the paper
  // adding 
  Point     xpq_hat  [3];    // three unit vector along x, y, z axes    
}; 


#endif // _M3D_H


