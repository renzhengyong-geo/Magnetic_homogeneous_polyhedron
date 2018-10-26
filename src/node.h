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

#ifndef _NODE_H 
#define _NODE_H

// C++   includes
#include <iomanip>
#include <string>
#include "point.h"

#define INVALID_UNIT 
class Node : public Point
{  
 public:
  Node (const double x=.0,
	const double y=.0,
	const double z=.0,
	const unsigned int id=0);
  Node (const Node& n);
  Node (const Point p,
	const unsigned int id=0);
  ~Node ();
  
  friend std::ostream& operator <<(std::ostream& os, const Node& node);
  
  // Overload operators.
  Node& operator=  (const Node& node);
  bool operator == (const Node& node);
  bool operator < (const Node& v) const; 
  
  void set_id(const unsigned int id) { _id=id; }
  unsigned int get_id() { return _id; }
 public:
  unsigned int               _id;
};


//-----------------------------------------------------------------------
inline
Node::Node (const double x, const double y,const double z, 
	    const unsigned int id) :
  Point    (x,y,z),
     _id    (id)
{ 
}

inline
Node::Node (const Node& n):
     Point (n),
     _id   (n._id)
{
}

inline
Node::Node (const Point p, 
	    const unsigned int id ):
  Point(p),
     _id  (id)
{  
}

inline
Node::~Node ()
{
}



inline
Node & Node::operator= (const Node& node)
{
  (*this)(0) = node(0);
  (*this)(1) = node(1);
  (*this)(2) = node(2); 
  _id        = node._id;

  return *this;
}

inline
bool Node::operator== (const Node& node) 
{
  if(_id!=node._id)
    return false;
  else 
    {
      //compare the coordinates.
      return Point(*this)==Point(node);
    }
}


inline
bool Node::operator < (const Node& rhs) const
{
  // Only compare xyz location.
  return (Point)(*this)<(Point)(rhs);
}



inline
std::ostream& operator <<(std::ostream& os, const Node& node)
{
  os<<node._id
    <<"\t"<<node._coords[0]
    <<"\t"<<node._coords[1]
    <<"\t"<<node._coords[2];
  return os;
}

#endif // NODE_H
