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

#include "grid.h"

Grid::Grid (std::string node_file, std::string tet_file):
  _node_file (node_file),
  _tet_file  (tet_file)
{
  // do the job
  this->read_node();
  this->read_tet();
}


Grid::~Grid ()
{
  // delete tets
  for(unsigned int i=0;i<(_elements.size());++i)
    if(_elements[i]!=NULL) {
      delete _elements[i]; _elements[i]=NULL;
    }
  this->_elements.clear();
  
  // delete the nodes.
  for(unsigned int i=0;i<(_nodes.size());++i) {
    if(_nodes[i]!=NULL) {
      delete _nodes[i];    _nodes[i]=NULL;
    }
  }
  this->_nodes.clear();  
  return;
}

/*
 get the i-th element
 */
Tet& Grid::get_elem (const unsigned int i) 
{
  assert (i < this->n_elems());
  assert (_elements[i] != NULL);
  return *this->_elements[i];
}


/*
 read nodes or sites
 */
void Grid::read_node()
{
  // Check input buffer
  std::ifstream node_stream(this->_node_file.c_str());
  assert(node_stream.good());
  
  unsigned int n_nodes=0;
  node_stream >> n_nodes; 
  assert(n_nodes>0);     
  this->_nodes.clear();
  this->_nodes.resize(n_nodes);
  
  for (unsigned int i=0; i<n_nodes; i++) {
    unsigned int node_lab=0;
    Point xyz;
    node_stream >> node_lab 
		>> xyz(0)    // x-coordinate value
		>> xyz(1)    // y-coordinate value
		>> xyz(2);   // z-coordinate value    
   Node* newnode = new Node(xyz, node_lab);
   this->_nodes[i]= newnode;
  }
  
  return;
}


/*
 read elements and magnetization vector M
 */
void Grid::read_tet ()
{
  // Check input buffer
  std::ifstream ele_stream(this->_tet_file.c_str());
  assert (ele_stream.good());
  assert (this->_nodes.size()>0);

  unsigned int n_elems=0;
  ele_stream >> n_elems;
  assert(n_elems>0);
  this->_elements.clear();
  this->_elements.resize(n_elems);
  
  for(unsigned int i=0; i<n_elems; i++) {
    unsigned int element_id;
    Node* v[4];
    ele_stream >> element_id;
    for (unsigned int j=0; j<4; j++)	{
      unsigned int node_id;
      ele_stream >> node_id;
      v[j] = this->_nodes[node_id];
    }
    double mx, my, mz;
    ele_stream >> mx >> my >> mz;
    Tet* tet = new Tet(element_id, v, Point(mx,my,mz));
    this->_elements[i] = tet;
  }
  
  this->_number_elems = n_elems; 
  
  return;
}




std::ostream& operator <<(std::ostream& os, Grid& mesh)
{
  os << "n_nodes()="    << mesh.n_nodes()
     << " \nn_elems()=" << mesh.n_elems()
     << " \n"; 
  return os;
}


/*
 Write the tetrahedral elements into VTk format so that 
 the grid can be shown by Paraview softeware.
 
 */
void Grid::write_out_vtk(const std::string name)
{ 
  // Open the of stream
  std::ofstream vtk_mesh((name+".vtk").c_str()); 
  if(!vtk_mesh.good()) {
    std::cerr<<"Can not open file:\t"<<name+".vtk"
	     <<std::endl;
  } else {
    vtk_mesh<<"# vtk DataFile Version 3.0\n" //File version and identifier
	    <<"Zhengyong Ren\n"
	    <<"ASCII\n";
    vtk_mesh<<"DATASET UNSTRUCTURED_GRID\n"; 
    vtk_mesh<<"\nPOINTS\t"<<this->_nodes.size()<<"\tdouble\n"; 
    for(unsigned int i=0; i<this->_nodes.size(); i++) {
      const Node& temp_node= (*this->_nodes[i]);
      vtk_mesh<<temp_node(0)<<"\t" 
	      <<temp_node(1)<<"\t"  
	      <<temp_node(2)<<"\n"; 
    }
    vtk_mesh<<"\nCELLS\t"<<this->n_elems()<<"\t"<<this->n_elems()*5<<"\n";
    for(unsigned int i=0; i<this->n_elems(); i++) {
      Tet& temp_elem= this->get_elem(i);
      vtk_mesh<<(unsigned int)4<<"\t" 
	      <<temp_elem.get_node_id(0)<<"\t"
	      <<temp_elem.get_node_id(1)<<"\t"
	      <<temp_elem.get_node_id(2)<<"\t" 
	      <<temp_elem.get_node_id(3)<<"\n";
    }
    vtk_mesh<<"\nCELL_TYPES\t"<<this->n_elems()<<"\n";
    for(unsigned int i=0; i<this->n_elems(); i++) {
      vtk_mesh<<(unsigned int)10<<"\n"; // 10--tet
    }
    vtk_mesh<<"\n";
               
  }
  
  vtk_mesh.close();
  // done
}


