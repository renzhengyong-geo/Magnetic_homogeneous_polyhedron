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

#ifndef  _GRID_H
#define  _GRID_H


#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>
#include "node.h"
#include "tet.h"

class Grid
{  
 public:
  Grid(std::string node_file, std::string tet_file);
  ~Grid();
  
  // tets
  unsigned int n_nodes() { return _nodes.size(); }
  unsigned int n_elems() { return _number_elems; }
  Tet&  get_elem (unsigned int i);
  
  void read_node();
  void read_tet ();
  
  friend std::ostream& operator <<(std::ostream& os, const Grid& mesh);
  void write_out_vtk(const std::string name);
  
 public:
  std::string                          _node_file;
  std::string                          _tet_file;
  std::vector<Node*>                   _nodes;
  std::vector<Tet*>	                   _elements;
  unsigned int                         _number_elems;  
};    

#endif // _GRID_H
