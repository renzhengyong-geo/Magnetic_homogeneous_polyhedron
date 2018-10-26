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


/*
	      0
   TET        
             /|\
          0 / | \
           /  |  \2
          /  1|   \
       1  ..4.|....3
          \   |   /
           \  |  /
          3 \ | /5
             \|/
              
              2
*/


#ifndef   _TET_H
#define   _TET_H

#include <cstdlib>
#include "node.h"

class Tet
{
public:
  Tet (unsigned int id, Node* v[4], Point M);
  unsigned int get_id ()    { return _id;}
  Point        get_M()      { return _M;}
  unsigned int get_node_id(unsigned int i);
  void         get_nodes(Node* nodes[4]);  
  void         get_facet(Node* facet[3],  unsigned int i);
  void         get_facet_normal_vector(Point& n,unsigned int i);
  void         get_facet_edge_nodes(Node* v[2], unsigned int i, unsigned int j);
  void         get_facet_edge_normal(Point& m, unsigned int i, unsigned int j);
  void         get_facet_edge_tangential(Point& t, unsigned int i, unsigned int j);
  void         get_facet_edge_ends (Node* a[2], unsigned int i, unsigned int j);
  double       get_size();  
 private:
  unsigned int       _id;
  Node*              _v[4];
  Point              _M;
};




inline
Tet::Tet(unsigned int id,
         Node*        v[4],
	 Point        M)
{
  this->_id = id;
  for(unsigned int i=0; i<4; i++) this->_v[i]=v[i];
  this->_M = M;
  
}

inline double Tet::get_size() 
{
  const Point & a=(*_v[1])+((*_v[0])*-1.0);
  const Point & b=(*_v[1])+((*_v[2])*-1.0);
  const Point & c=(*_v[1])+((*_v[3])*-1.0);
  return std::abs( a*(c.cross(b)) )/6.;
}

inline
unsigned int Tet::get_node_id(unsigned int i)
{
  assert(i<4);
  return this->_v[i]->get_id();
}


inline  
void Tet::get_nodes(Node* nodes[4])
{
  for(int i=0; i<4; i++) nodes[i] = this->_v[i];
  return;
}

inline
void Tet::get_facet(Node* facet[3],  unsigned int i)
{
  assert(i<4);
  switch(i) {
  case 0:
    facet[0] = _v[1]; facet[1] = _v[2]; facet[2] = _v[3]; 
    break;
  case 1:
    facet[0] = _v[0]; facet[1] = _v[2]; facet[2] = _v[3]; 
    break;
  case 2:
    facet[0] = _v[0]; facet[1] = _v[1]; facet[2] = _v[3]; 
    break;
  case 3:
    facet[0] = _v[0]; facet[1] = _v[1]; facet[2] = _v[2];   
    break;
  default:
    std::abort();
  }
  return;
}


inline
void Tet::get_facet_normal_vector(Point& n,unsigned int i)
{
  Point normal, facet_center;
  assert(i<4);
  switch(i) {
   case 0:
     normal    = ((*_v[3])-(*_v[2])).cross((*_v[1])-(*_v[2])).unit();
     facet_center = ((*_v[1])+(*_v[2])+(*_v[3]))*(1.0/3.0);
     break;
  case 1:
    normal    = ((*_v[0])-(*_v[3])).cross((*_v[2])-(*_v[3])).unit();
    facet_center = ((*_v[0])+(*_v[2])+(*_v[3]))*(1.0/3.0);
    break;
  case 2:
    normal = ((*_v[3])-(*_v[0])).cross((*_v[1])-(*_v[0])).unit();
    facet_center = ((*_v[0])+(*_v[1])+(*_v[3]))*(1.0/3.0);
    break;
  case 3:
    normal = ((*_v[2])-(*_v[1])).cross((*_v[0])-(*_v[1])).unit();
    facet_center = ((*_v[0])+(*_v[1])+(*_v[2]))*(1.0/3.0);
    break;
  default:
    std::abort();
  }
  
  Point center = ((*_v[0])+(*_v[1])+(*_v[2])+(*_v[3]))*(0.25);
  // we need outgoing normal vector
  if(normal*(facet_center-center)<0.0) normal= normal*-1.0;
  
  n = normal;
  return;
  
}


/*           1
        0---------2
          \      /
           \    /
          2 \  / 0
             \/              
              1
*/
inline
void  Tet::get_facet_edge_nodes(Node* v[2], unsigned int i, unsigned int j)
{
  assert(i<4);
  assert(j<3);
  Node* facet[3];
  this->get_facet(facet, i);
  switch(j) {
  case 0:
    v[0] = facet[1];    v[1] = facet[2];
    break;
  case 1:
    v[0] = facet[2];    v[1] = facet[0];
    break;
  case 2:
    v[0] = facet[0];    v[1] = facet[1];
    break;
  default:
    std::abort();
  }
  return;
}


inline
void  Tet::get_facet_edge_normal(Point& m, unsigned int i, unsigned int j)
{
  assert(i<4);
  assert(j<3);
  Node* v[2];
  this->get_facet_edge_nodes(v, i, j);
  Point t = ((*v[1])-(*v[0])).unit();
  Point n;
  this->get_facet_normal_vector(n, i);
  Point temp_m = t.cross(n);
  // we need outgoing normal on edge
  Node* facet[3];
  this->get_facet(facet,  i);
  Point facet_center =  ((*facet[0])+(*facet[1])+(*facet[2]))*(1.0/3.0);
  Point edge_center   =  ((*v[1])+(*v[0]))*0.5;
  if(temp_m*(edge_center-facet_center)<0.0) temp_m= temp_m*-1.0;
   
  m = temp_m;
  return;
}

inline
void Tet::get_facet_edge_tangential(Point& t, unsigned int i, unsigned int j)
{
  assert(i<4);
  assert(j<3);
  Point m;
  this->get_facet_edge_normal(m, i, j);
  Point n;
  this->get_facet_normal_vector(n, i);
  t = n.cross(m);
}



/*
a[0]-----------a[1]
(a-)           (a+)
*/
inline
void  Tet::get_facet_edge_ends (Node* a[2], unsigned int i, unsigned int j)
{
  assert(i<4);
  assert(j<3);
  Point t;
  this->get_facet_edge_tangential(t, i, j);
  Node* v[2];
  this->get_facet_edge_nodes(v, i, j);

  Point temp_t = ((*v[1])-(*v[0])).unit();
  if(t*temp_t>0.0) {
    a[0] = v[0];
    a[1] = v[1];
  }else if(t*temp_t<0.) {
    a[0] = v[1];
    a[1] = v[0];
  }
  return;
}



#endif     // _Tet_h


