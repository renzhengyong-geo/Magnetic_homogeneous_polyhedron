/*

author: Chen Chaojian, ETH Zurich, 2017

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>


/*

Heath, P. J., S. Greenhalgh, and N. G. Direen, 2005, Modelling
gravity and magnetic gradient tensor responses for explo-
ration within the regolith: Exploration Geophysics, 36, 357–
364.
*/



double factor[8] = {1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0};
int index[8][3]  = {{2,2,2},{1,2,2},{2,1,2},{1,1,2},
		    {2,2,1},{1,2,1},{2,1,1},{1,1,1}};



// Bxx
double Bxx_kernel(double x, double y, double z)
{
  double r    = std::sqrt(x*x+y*y+z*z);
  double temp = x/(r*(y+r));
  return temp;
}


double Bxx(double x1, double x2,
	   double y1, double y2,
	   double z1, double z2,
	   double Fe, double k)
{
  double x[2]={x1,x2};
  double y[2]={y1,y2};
  double z[2]={z1,z2};
  
  // double S   = (Fe*k)/(6.673e-11*m);
  double S   = (Fe*k);
  double sum=0.;
  for(int i=0; i<8; i++) {
    int index_x = index[i][0]-1;
    int index_y = index[i][1]-1;
    int index_z = index[i][2]-1;
    sum+=  Bxx_kernel(x[index_x],y[index_y],z[index_z])*factor[i];
  }
  return sum*S;
}

// Bxy
double Bxy_kernel(double x, double y, double z)
{
  double r    = std::sqrt(x*x+y*y+z*z);
  double temp = 1./r;
  return temp;
}


double Bxy(double x1, double x2,
	   double y1, double y2,
	   double z1, double z2,
	   double Fe, double k)
{
  double x[2]={x1,x2};
  double y[2]={y1,y2};
  double z[2]={z1,z2};
  
  // double S   = (Fe*k)/(6.673e-11*m);
  double S   = (Fe*k);
  double sum=0.;
  for(int i=0; i<8; i++) {
    int index_x = index[i][0]-1;
    int index_y = index[i][1]-1;
    int index_z = index[i][2]-1;
    sum+=  Bxy_kernel(x[index_x],y[index_y],z[index_z])*factor[i];
  }
  return sum*S;
}

// Bxz
double Bxz_kernel(double x, double y, double z)
{
  double r    = std::sqrt(x*x+y*y+z*z);
  double temp = z/(r*(y+r));
  return temp;
}


double Bxz(double x1, double x2,
	   double y1, double y2,
	   double z1, double z2,
	   double Fe, double k)
{
  double x[2]={x1,x2};
  double y[2]={y1,y2};
  double z[2]={z1,z2};
  
  // double S   = (Fe*k)/(6.673e-11*m);
  double S   = (Fe*k);
  double sum=0.;
  for(int i=0; i<8; i++) {
    int index_x = index[i][0]-1;
    int index_y = index[i][1]-1;
    int index_z = index[i][2]-1;
    sum+=  Bxz_kernel(x[index_x],y[index_y],z[index_z])*factor[i];
  }
  return sum*S;
}

// Byy
double Byy_kernel(double x, double y, double z)
{
  double r    = std::sqrt(x*x+y*y+z*z);
  double temp = y/(r*(x+r));
  return temp;
}


double Byy(double x1, double x2,
	   double y1, double y2,
	   double z1, double z2,
	   double Fe, double k)
{
  double x[2]={x1,x2};
  double y[2]={y1,y2};
  double z[2]={z1,z2};
  
  // double S   = (Fe*k)/(6.673e-11*m);
  double S   = (Fe*k);
  double sum=0.;
  for(int i=0; i<8; i++) {
    int index_x = index[i][0]-1;
    int index_y = index[i][1]-1;
    int index_z = index[i][2]-1;
    sum+=  Byy_kernel(x[index_x],y[index_y],z[index_z])*factor[i];
  }
  return sum*S;
}


// Byz
double Byz_kernel(double x, double y, double z)
{
  double r    = std::sqrt(x*x+y*y+z*z);
  double temp = z/(r*(x+r));
  return temp;
}


double Byz(double x1, double x2,
	   double y1, double y2,
	   double z1, double z2,
	   double Fe, double k)
{
  double x[2]={x1,x2};
  double y[2]={y1,y2};
  double z[2]={z1,z2};
  
  // double S   = (Fe*k)/(6.673e-11*m);
  double S   = (Fe*k);
  double sum=0.;
  for(int i=0; i<8; i++) {
    int index_x = index[i][0]-1;
    int index_y = index[i][1]-1;
    int index_z = index[i][2]-1;
    sum+=  Byz_kernel(x[index_x],y[index_y],z[index_z])*factor[i];
  }
  return sum*S;
}


int main() 
{
  double x[2] ={-5,5};
  double y[2] ={-5,5};
  double z[2] ={ 0,-10};
  double Fe = 1.;        // magnitude of the including magnetic field
  double k  = 200.;      // the magnetic susceptibility of the rectangular prism
  

  int n=21;
  std::vector<double> BXX(n), BXY(n), BXZ(n), BYY(n), BYZ(n), BZZ(n);
  std::ofstream data("result.dat");
  for(int i=0; i<n; i++) {
    // observation point 
    // (x0, y0, z0)
    double x0 = 6.;
    double y0 = -25. + 2.5*i;
    double z0 = 0.;
    // move original point onto point(x0, y0, z0)
    double X1 = x[0]- x0;
    double X2 = x[1]- x0;
    double Y1 = y[0]- y0;
    double Y2 = y[1]- y0;
    double Z1 = z[0]- z0;
    double Z2 = z[1]- z0;
	// Unit: nT/s
    BXX[i] = 100*Bxx(X1, X2, Y1, Y2, Z1, Z2, Fe, k);
    BXY[i] = 100*Bxy(X1, X2, Y1, Y2, Z1, Z2, Fe, k);
    BXZ[i] = 100*Bxz(X1, X2, Y1, Y2, Z1, Z2, Fe, k);
    BYY[i] = 100*Byy(X1, X2, Y1, Y2, Z1, Z2, Fe, k);
    BYZ[i] = 100*Byz(X1, X2, Y1, Y2, Z1, Z2, Fe, k);
    BZZ[i] = -1.*(BXX[i] + BYY[i]);

	data<<std::setprecision(20)<<BXX[i]<<"  "
            <<std::setprecision(20)<<BXY[i]<<"  "
            <<std::setprecision(20)<<BXZ[i]<<"  "
		<<std::setprecision(20)<<BXY[i]<<"  "
            <<std::setprecision(20)<<BYY[i]<<"  "
            <<std::setprecision(20)<<BYZ[i]<<"  "
		<<std::setprecision(20)<<BXZ[i]<<"  "
		<<std::setprecision(20)<<BYZ[i]<<"  "
            <<std::setprecision(20)<<BZZ[i]<<std::endl; 

  }

  data.close();

  return 0;
}
