#include "Eigen/Core"
#include <iostream>
#include <vector>
#include "Eigen/src/Core/Matrix.h"
#include "calc_factor_group.cpp"
//#include "crystal_class.cpp"
#include "io_func.cpp"

#define PREC 1e-6

int main(int argc, char *argv[])
{
/*reads in lattice and basis from a POSCAR
 * find point group of lattice and then factor group of crystal
 *prints factor group to screen.
 */
  if (argc !=2){
      std::cout<<"WRONG NUMER OF INPUT ARGUMENTS!"<<std::endl;
      std::cout<<"Please specify a path to a POSCAR file afer executable."<<std::endl;
      return 1;
  }
  Crystal_Structure my_struct = read_poscar(argv[1]);
//  Eigen::Matrix3d  lattice;
//  std::vector<Eigen::Vector3d> basis;
//  Crystal_Structure my_struct(lattice, basis);/  
//  my_struct.lattice.col(0) <<1.0, 1.0, -1.0;
//  my_struct.lattice.col(1) <<-1.0, 1.0, 1.0;
//  my_struct.lattice.col(2) << 1.0, -1.0, 1.0;
//  Eigen::Vector3d coord; coord <<0.0, 0.0, 0.0;
//  my_struct.basis.push_back(coord);
 // coord << 0.5, 0.5, 0.5;
//  my_struct.basis.push_back(coord);
//  coord << 0.5, 0.0, 0.5;
//  my_struct.basis.push_back(coord);
//  coord << 0.5, 0.5, 0.0;
//  my_struct.basis.push_back(coord);
  int len_basis = my_struct.basis.size();
  std::cout<<"This is the lattice: \n"<<my_struct.lattice<<std::endl;
  std::cout<<"This is the basis: \n";
  for(int i=0; i<len_basis; i++){
      my_struct.basis[i]= fractional_to_cartesian(my_struct.lattice, my_struct.basis[i]);
      std::cout<<my_struct.basis[i]<<std::endl;
  }


  std::cout<<"MY Struct lattice is : \n" << my_struct.lattice<<std::endl;
  auto pt_group = calc_point_group(my_struct.lattice);
  print_ptgroup(pt_group, my_struct.lattice);
  auto my_factor_group = find_factor_group(my_struct);
  print_factor_group(my_factor_group);
  


  return 0;
}

