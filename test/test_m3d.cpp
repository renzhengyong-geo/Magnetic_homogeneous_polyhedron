/* m3d
 * 	These files may be found at:
 * 	http://software.seg.org/2017/0006
 *
 * 	Please first read and accept:
 * 	http://software.seg.org/disclaimer.txt
*/

/*
  Zhengyong Ren, zhengyong@csu.edu.cn
  Central South University of China, 2018
*/


#define NDEBUG
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <fstream>
#include "grid.h"
#include "m3d.h"


//------------------read----------------------------------------------
void read_sites(std::string site_file, std::vector<Point>& sites)
{
  std::ifstream site_stream(site_file.c_str());  
  unsigned int n_sites=0;
  site_stream >> n_sites; 
  assert(n_sites>0);     
  sites.clear();
  sites.resize(n_sites);  
  for (unsigned int i=0; i<n_sites; i++) {
    unsigned int site_id=0;
    Point xyz;
    site_stream >> site_id 
		>> xyz(0)    // x-coordinate value
		>> xyz(1)    // y-coordinate value
		>> xyz(2);   // z-coordinate value    
    sites[i]= xyz;
  }
  return;
}

//------------------write----------------------------------------------
void write_vector(std::string name, std::vector<double>& v)
{
  // Open the of stream
  std::ofstream out(name.c_str()); 
  if(!out.good()) {
    std::cerr<<"Can not open file:\t"<<name
	     <<std::endl;
  } else {
    for(unsigned int i=0; i<v.size(); i++)
      out<<std::setprecision(16)<<v[i]<<"\n";
  }
  out.close();    
  
}

void write_vector(std::string name, std::vector<Point>& v)
{
  // Open the of stream
  std::ofstream out(name.c_str()); 
  if(!out.good()) {
    std::cerr<<"Can not open file:\t"<<name
	     <<std::endl;
  } else {
    for(unsigned int i=0; i<v.size(); i++)
      out<<std::setprecision(16)<<v[i](0)<<"\t"<<v[i](1)<<"\t"<<v[i](2)<<"\n";
  }
  out.close();    
  
}

void write_vector(std::string name, std::vector<Point> & TBx, 
		       std::vector<Point> & TBy, 
		       std::vector<Point> & TBz)
{
  // Open the of stream
  std::ofstream out(name.c_str()); 
  if(!out.good()) {
    std::cerr<<"Can not open file:\t"<<name
	     <<std::endl;
  } else {
    for(unsigned int i=0; i<TBx.size(); i++) {
      out<<std::setprecision(16)<<TBx[i](0)<<"\t"<<TBx[i](1)<<"\t"<<TBx[i](2)<<"\t";
      out<<std::setprecision(16)<<TBy[i](0)<<"\t"<<TBy[i](1)<<"\t"<<TBy[i](2)<<"\t";
      out<<std::setprecision(16)<<TBz[i](0)<<"\t"<<TBz[i](1)<<"\t"<<TBz[i](2)<<"\n";
    }
  }
  out.close(); 
  return;
}


//---------------------------------------------------------------------------
int main(int argc, char *argv[]) 
{
  if (argc != 4) {
    printf("Usage: %s node_file tet_file site_fie\n",argv[0]);
    return 1;    
  }  
  std::string node_file(argv[1]);
  std::string tet_file (argv[2]);
  std::string site_file(argv[3]);

  // read grid
  Grid grid(node_file, tet_file);
  // read sites
  std::vector<Point> sites;
  read_sites(site_file, sites);
  write_vector("xyz.dat", sites);
    
  // loop each site and tet pair
  const int n_sites = sites.size();
  std::vector<double> Vs;
  std::vector<Point>  Bs;  
  std::vector<Point> Tx, Ty, Tz;
  Vs.resize(n_sites);
  Bs.resize(n_sites);
  Tx.resize(n_sites);
  Ty.resize(n_sites);
  Tz.resize(n_sites);

  // closed-form solutions
#pragma omp parallel for
  for(int i=0; i<n_sites; i++) {
    Point& r = sites[i];
    double temp_v; Point temp_b(0,0,0); double temp_tb[3][3];
    temp_v=0;
    for(int l=0; l<3; l++)
      for(int k=0; k<3; k++)
        temp_tb[l][k]=0.;
    
    for(int j=0; j<grid.n_elems(); j++) {
      Tet& tet = grid.get_elem(j);
      M3D m3d(tet, r);
      double v; Point b; double tb[3][3];
      m3d.V(v);
      m3d.B(b);
      m3d.T(tb);
      temp_v += v;
      temp_b = temp_b + b;
      for(int l=0; l<3; l++)
        for(int k=0; k<3; k++)
         temp_tb[l][k] += tb[l][k];
    }
    Vs[i] = temp_v;
    Bs[i] = temp_b*1e9;
    Tx[i]= Point(temp_tb[0][0], temp_tb[0][1], temp_tb[0][2])*1e9;
    Ty[i]= Point(temp_tb[1][0], temp_tb[1][1], temp_tb[1][2])*1e9;
    Tz[i]= Point(temp_tb[2][0], temp_tb[2][1], temp_tb[2][2])*1e9;
  } // done. 


  write_vector("V.dat",   Vs);         // NA-1m-3
  write_vector("B.dat",   Bs);         // nT
  write_vector("T.dat",  Tx, Ty, Tz);  // nT/m

  return 0;  
}






