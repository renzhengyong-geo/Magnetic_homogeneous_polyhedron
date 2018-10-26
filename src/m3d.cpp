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


/*
 Input: a tetrahedon and an observation site
 Output: ni_hat  (outgoing normal vector of 4 surfaces)
         hi (heights from site to 4 surfaces
         Rij0 and Rij1
         oi (projection point of site on 4 surfaces)
         Sij0 and Sij1
         Rij (distance from site to edges of 4 surfaces)
         eij_hat (tangential normal vector on edges)
         mij_hat (outgoing normal vector of edges)
 */


// C++ includes
#include "m3d.h"
#include <iostream>

M3D::M3D(Tet& tet_, Point r_):
  tet   (tet_),
  r     (r_),
  TOL   (1e-14),
  Cm    (1e-07)  
{
  // calculate variables for facets and edges
  for(int i=0; i<4; i++) {
     Point n;
     this->tet.get_facet_normal_vector(n, i);
     ni_hat[i] = n;
     double w0; 
     for(int j=0; j<3; j++) {
        Node* a[2];
        this->tet.get_facet_edge_ends (a, i, j);
        Rij0[i][j] = (this->r-(*a[0])).size();
        Rij1[i][j]  = (this->r-(*a[1])).size();
        if(j==0) { 
        w0 = (r-(*a[0]))*n;}
     }
     hi[i]=w0;
 
     Point o = r - n*w0;
     oi[i] = o;
     for(int j=0; j<3; j++) {
        Point m, e; Node* a[2];
        this->tet.get_facet_edge_normal(m, i, j);
        this->tet.get_facet_edge_tangential(e, i, j);
        this->tet.get_facet_edge_ends (a, i, j);
        Sij0[i][j] = ((*a[0])-o)*e;
        Sij1[i][j]  = ((*a[1])-o)*e;
        mij[i][j]= ((*a[0])-o)*m;
        Rij[i][j] = std::sqrt(w0*w0+mij[i][j]*mij[i][j]);
        mij_hat[i][j]  = m;
        eij_hat[i][j]  = e; //e=n.cross(m)
        // for added two edge variables:
        Point temp_rij_n[4][3];
        temp_rij_n[i][j] = o + mij_hat[i][j]*mij[i][j];
        if(std::abs(mij[i][j])<TOL)
          rhoij_n[i][j] = temp_rij_n[i][j]-temp_rij_n[i][j];// rhoij_n[i][j]=(0,0,0)
        else
          rhoij_n[i][j] = (temp_rij_n[i][j]-o).unit();
        rij_n[i][j] = temp_rij_n[i][j]-r;// in local coordinate with r_local=(0,0,0), which is assumed in the manuscript.
     }
  }

  Point x_hat(1.,0,0);
  Point y_hat(0,1.,0);
  Point z_hat(0,0,1.);
  xpq_hat[0] = x_hat;
  xpq_hat[1] = y_hat; 
  xpq_hat[2] = z_hat; 
   
}


/********************************************************************/
/**************calculating different kinds of integrals**************/
/********************************************************************/

/*
 Input: surface index i 
 Output: k0=\iint_{T_{i}} 1/R
 */
void M3D::compute_1_R(unsigned int i, double &k0)
{
  double w0 = hi[i];
  double temp_k=0.;
  
  for(int j=0; j<3; j++) {
    double part1=0, part2=0;
    if(std::abs(w0)>TOL) {
      double beta = std::atan((mij[i][j]*Sij1[i][j])/(Rij[i][j]*Rij[i][j]+std::abs(w0)*Rij1[i][j]))
	- std::atan((mij[i][j]*Sij0[i][j])/(Rij[i][j]*Rij[i][j]+std::abs(w0)*Rij0[i][j]));
      part2 = std::abs(w0)*beta;
    }
    double mj0 = mij[i][j];
    if(std::abs(mj0)>TOL) {
      part1 = mj0*std::log((Rij1[i][j]+Sij1[i][j])/(Rij0[i][j]+Sij0[i][j]));
    }
    temp_k += part1 - part2;
  }
  
  k0 = temp_k;
  
  return;
  
}
       
/*
 Input: surface index i
 Output: Ki=\iint_{T_{i}} \nabla 1/ R 
 where \nabla is on the observation site r
 */

void M3D::compute_grad_1_R(unsigned int i, Point &Ki)
{
  Point n   = ni_hat[i];
  double w0 = hi[i];
  Point temp_k_minus;
  
  for(int j=0; j<3; j++) {
    Point part1, part2;
    if(std::abs(w0)>TOL) {
      double betaij= atan((mij[i][j]*Sij1[i][j])/(Rij[i][j]*Rij[i][j]+std::abs(w0)*Rij1[i][j]))
	- atan((mij[i][j]*Sij0[i][j])/(Rij[i][j]*Rij[i][j]+std::abs(w0)*Rij0[i][j]));
      part2  = n*sgn(w0)*betaij;
    }
    double Aij=0.;
    if(std::abs(Rij[i][j])>TOL) {
      Aij = std::log((long double)(Rij1[i][j]+Sij1[i][j])) - std::log((long double)(Rij0[i][j]+Sij0[i][j]));
    }else if(Sij1[i][j]>0.&&Sij0[i][j]>0.) {
      Aij = std::log(Sij1[i][j]) - std::log(Sij0[i][j]);
    }else if(Sij1[i][j]<0.&&Sij0[i][j]<0.) {
      Aij = std::log(Sij0[i][j]*-1.0)- std::log(Sij1[i][j]*-1.0);
    }else {
      std::cout<<"Site must not be on any edge!\n";
      std::abort();
    }
    part1 = mij_hat[i][j]*Aij;
    temp_k_minus = temp_k_minus + part1 + part2;
  }
  
  Ki = temp_k_minus*-1.0;
  
  return;
  
}

/*
 Input: 4 nodes in vector T
 Output: volume(size) of tetrahedron
 */
double M3D::get_tri_size(std::vector<Point>& T)
{
  assert(T.size()==3);
  Point a=T[1]+(T[0]*-1.0);
  Point b=T[2]+(T[0]*-1.0);
  Point c=a.cross(b);
  
  return c.size()*0.5;
}

/*
 Input: surface index i
 Output: k2=\iint_{T_{i}} \nabla \nabla 1/ R
 where \nabla is on the observation site r
 */
void M3D::compute_grad_grad_1_R(unsigned int i, double k2[3][3])
{
  double temp_k2[3][3];
  this->set(temp_k2, 0.);

  for(int p=0; p<3; p++)
      for(int q=0; q<3; q++)
        {
          double kxpxqR5 = 0.;
          this->compute_xpxq_R5(i,p,q,kxpxqR5);
          if(p==q)
          { 
            double k1R3=0.;
            this->compute_1_R3(i,k1R3);
            temp_k2[p][q] = 3.*kxpxqR5 - k1R3;
          }
          else
          {
            temp_k2[p][q] = 3.*kxpxqR5;
          }
        } 

  // assign k2
  for(int j=0; j<3; j++)
    for(int k=0; k<3; k++) 
      k2[j][k] = temp_k2[j][k];

  return;
}

/*
 Input: surface index i
 Output: k3=\iint_{1/R^3} 
 */

void M3D::compute_1_R3(unsigned int i, double &k3)
{
  double temp_k3=0.;
  for(int j=0; j<3; j++)
  {
     // mk1Rr2=mij*\int_{1/(R*\_rho^2)}
     double lmijr2R=0.;
     this->line_mij_r2R(i,j,lmijr2R);
     temp_k3 -= lmijr2R;
  }

  double belta=0.; // the angle on the polygon of oi
  for(int k=0; k<3; k++)
  {
    double arc1,arc0=0.;
    if(std::abs(mij[i][k])<TOL)
      belta += 0.; //betta_ij = 0
    else
    {
      arc1 = std::atan(Sij1[i][k]/std::abs(mij[i][k]));
      arc0 = std::atan(Sij0[i][k]/std::abs(mij[i][k]));   
      belta += (mij_hat[i][k]*rhoij_n[i][k])*(arc1-arc0);
    }
  }

  if(std::abs(hi[i])<TOL)
    temp_k3 += 0.;
  else
    temp_k3 += belta/std::abs(hi[i]);

  // assign k3
  k3 = temp_k3;
  return;
}

/*
 Input: surface index i, xp index p and xq index q
 Output: k4=\iint_{xp*xq/R^5} 
 where \nabla is on the observation site r
 */
void M3D::compute_xpxq_R5(unsigned int i, unsigned int p, unsigned int q, double &k4)
{
  double temp_k4 = 0.;
  double k3,k5 = 0.;

  for(int j=0; j<3; j++)
  {
    double l1R3,lxqR3=0.; //lxqR3=\int_{xq/R^3}
    this->line_1_R3(i,j,l1R3);
    lxqR3 = (rij_n[i][j]*xpq_hat[q])*l1R3;
    lxqR3 += (eij_hat[i][j]*xpq_hat[q])*(1./Rij0[i][j]-1./Rij1[i][j]);   
    temp_k4 -= 1./3.*(mij_hat[i][j]*xpq_hat[p])*lxqR3;
  }
  this->compute_xq_R5(i, q, k5);
  temp_k4 -= (xpq_hat[p]*ni_hat[i])*hi[i]*k5;
  this->compute_1_R3(i, k3);
  temp_k4 += ((xpq_hat[p]*xpq_hat[q]-(xpq_hat[p]*ni_hat[i])*(ni_hat[i]*xpq_hat[q]))/3.)*k3;
  // assign k4
  k4 = temp_k4;
  return;  
}

/*
 Input: surface index i, xq index q
 Output: k5=\iint_{xq/R^5} 
 */
void M3D::compute_xq_R5(unsigned int i, unsigned int q, double &k5)
{
  double temp_k5,k6=0.;
  for(int j=0; j<3; j++)
  {
    double l1R3=0.;
    this->line_1_R3(i,j,l1R3);
    temp_k5 -= 1./3.*(mij_hat[i][j]*xpq_hat[q])*l1R3;
      }
  this->compute_1_R5(i, k6);
  temp_k5 -= (xpq_hat[q]*ni_hat[i])*hi[i]*k6;
  
  // assign k5
  k5 = temp_k5;
  return;
}

/*
 Input: surface index i
 Output: k6=\iint_{1/R^5} 
 */
void M3D::compute_1_R5(unsigned int i, double &k6)
{
  double temp_k6=0.;
  for(int j=0; j<3; j++)
  {
     double lmijr2R3=0.;
     this->line_mij_r2R3(i,j,lmijr2R3);     
     temp_k6 -= lmijr2R3/3.;
  }
  
  double belta=0.; // the angle on the polygon of oi
  for(int k=0; k<3; k++)
  {
    double arc1,arc0=0.;
    if(std::abs(mij[i][k])<TOL)
      belta += 0.; //betta_ij = 0
    else
    {
      arc1 = std::atan(Sij1[i][k]/std::abs(mij[i][k]));
      arc0 = std::atan(Sij0[i][k]/std::abs(mij[i][k])); 
      belta += (mij_hat[i][k]*rhoij_n[i][k])*(arc1-arc0);
    }
  }

  double temp_h=0.;
  temp_h = std::abs(hi[i]);
  if(temp_h<TOL)
    temp_k6 += 0.;
  else
  temp_k6 += belta/(3.*temp_h*temp_h*temp_h);

  // assign k6
  k6 = temp_k6;
  return;
}

/*
 Input: surface index i, line index j
 Output: k7=\int_{1/R^3} 
 */
void M3D::line_1_R3(unsigned int i, unsigned int j, double &k7)
{
  double temp_k7=0.;
  if((std::abs(hi[i])>TOL)||(std::abs(mij[i][j])>TOL)) // hi !=0 or mij !=0
  {
    temp_k7 = Sij1[i][j]/((hi[i]*hi[i]+mij[i][j]*mij[i][j])*Rij1[i][j]);
    temp_k7 -= Sij0[i][j]/((hi[i]*hi[i]+mij[i][j]*mij[i][j])*Rij0[i][j]);
  }
  else if((Sij0[i][j]*Sij1[i][j])>0.) // hi =0, mij = 0, r is not on the edg Cij.
  {
    double sub1,sub0=0.;
    sub1 = 1./(2.*Sij1[i][j]*Sij1[i][j]);
    sub0 = 1./(2.*Sij0[i][j]*Sij0[i][j]);
    temp_k7 = std::abs(sub1-sub0);
  }
  else
  {
    std::cout<<"Site must not be on any edge!\n";
    std::abort();
  }

  // assign k7
  k7 = temp_k7;
  return;
}

/*
 Input: surface index i, line index j
 Output: k8=\int_{m_{ij}/(rho^2*R)} 
 */
void M3D::line_mij_r2R(unsigned int i, unsigned int j, double &k8)
{
  double temp_k8=0.;
  if((std::abs(mij[i][j])>TOL)&&(std::abs(hi[i])>TOL)) //mij !=0 and hi !=0
  {
    double arc1,arc0=0.;
    arc1 = std::atan((std::abs(hi[i])*Sij1[i][j])/(mij[i][j]*Rij1[i][j]));
    arc0 = std::atan((std::abs(hi[i])*Sij0[i][j])/(mij[i][j]*Rij0[i][j]));
    temp_k8 = (arc1-arc0)/std::abs(hi[i]);
  }
  else if((std::abs(mij[i][j])>TOL)&&(std::abs(hi[i])<TOL)) //mij !=0 and hi=0 
  {
    double div1,div0=0.;
    div1 = Sij1[i][j]/std::sqrt(Sij1[i][j]*Sij1[i][j]+mij[i][j]*mij[i][j]);  
    div0 = Sij0[i][j]/std::sqrt(Sij0[i][j]*Sij0[i][j]+mij[i][j]*mij[i][j]); 
    temp_k8 = (div1-div0)/mij[i][j];
  }
  else //mij =0
  {
    temp_k8 = 0.;
  }

  //asign k8
  k8 = temp_k8;
  return;
}

/*
 Input: surface index i, line index j
 Output: k9=\int_{m_{ij}/(rho^2*R^3)} 
 */
void M3D::line_mij_r2R3(unsigned int i, unsigned int j, double &k9)
{
  double temp_k9=0.;
  if((std::abs(mij[i][j])>TOL)&&(std::abs(hi[i])>TOL)) //mij !=0 and hi !=0
  {
    double arc1,arc0=0.;
    double div1,div0=0.;
    arc1 = std::atan((std::abs(hi[i])*Sij1[i][j])/(mij[i][j]*Rij1[i][j]));
    arc0 = std::atan((std::abs(hi[i])*Sij0[i][j])/(mij[i][j]*Rij0[i][j]));
    div1 = Sij1[i][j]/Rij1[i][j];
    div0 = Sij0[i][j]/Rij0[i][j];
    temp_k9 = (arc1-arc0)/(std::abs(hi[i])*std::abs(hi[i])*std::abs(hi[i]));
    temp_k9 -= mij[i][j]*(div1-div0)/(hi[i]*hi[i]*(hi[i]*hi[i]+mij[i][j]*mij[i][j])); 
  }
  else if((std::abs(mij[i][j])>TOL)&&(std::abs(hi[i])<TOL)) //mij !=0 and hi=0 
  {
    double sqr1,sqr0,sub1,sub2=0.;
    sqr1 = std::sqrt(Sij1[i][j]*Sij1[i][j]+mij[i][j]*mij[i][j]);
    sqr0 = std::sqrt(Sij0[i][j]*Sij0[i][j]+mij[i][j]*mij[i][j]);
    sub1 = Sij1[i][j]/(sqr1*sqr1*sqr1)-Sij0[i][j]/(sqr0*sqr0*sqr0);
    sub2 = Sij1[i][j]/sqr1-Sij0[i][j]/sqr0;
    temp_k9 = sub1/(3.*mij[i][j])+2.*sub2/(3.*mij[i][j]*mij[i][j]*mij[i][j]);
  }
  else //mij =0
  {
    temp_k9 = 0.;
  }

  // assign k9
  k9 = temp_k9; 
  return;
}

// base function used in surface integral computation
/*
 a=b
 */
void M3D::set(double a[3][3], double b)
{
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++) 
      a[i][j]=b;
  return;
}


/*
 sign of real number m
 */
double M3D::sgn(const double m) 
{
  if (std::abs(m)<TOL) return 0.; 
  else return m/std::abs(m);
}


/********************************************************************/
/**************calculating Magentic anomalies***************/
/********************************************************************/

/*
 calculate the magnetic potential V
 */
void M3D::V(double& v)
{
  double temp_V=0.;
  Point M  = this->tet.get_M();

  for(int i=0; i<4; i++) {
    double k0;
    this->compute_1_R(i, k0); 
    double mdotn = M*ni_hat[i];
    temp_V+= k0*mdotn;
  }
  v = temp_V*Cm;  

  return;  
}


/*
 calculate the magnetic field B
 */
void M3D::B(Point& b)
{
  Point temp_B(0,0,0);
  Point M  = this->tet.get_M();

  for(int i=0; i<4; i++) {
    Point k1(0,0,0);
    this->compute_grad_1_R(i, k1);  
    double mdotn = M*ni_hat[i];
    temp_B = temp_B + k1*mdotn;
  }
  b = temp_B*(Cm*-1.0);  

  return ;   
}

/*
 calculate the magnetic gradient tensor T
 */
void M3D::T(double tb[3][3])
{
  double temp_tb[3][3];
  this->set(temp_tb, 0);
  Point M  = this->tet.get_M();

  // calculating tensor temp_tb
  for(int i=0; i<4; i++)
     {
        double k2[3][3];
        this->set(k2, 0);
        this-> compute_grad_grad_1_R(i,k2);
        double ndotM = ni_hat[i]*M;
        for(int j=0; j<3; j++)
             for(int k=0; k<3; k++)
                temp_tb[j][k] -= 1.0*Cm*ndotM*k2[j][k];         
     }

  // tb = temp_tb
  for(int i=0; i<3; i++)
     for(int j=0; j<3; j++) 
        tb[i][j] = temp_tb[i][j];

  return ;   
}

