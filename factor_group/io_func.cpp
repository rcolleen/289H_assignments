#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include "Eigen/src/Core/util/ForwardDeclarations.h"
#include "symop.cpp" //comment out to compile factor group
#include "crystal_class.cpp"//comment out to compile factor group 
#include "categorize_symop_main.cpp"

using namespace std;

auto read_lattice(string filename)
{
//read file and take each word and use it to populate L, an Eigen matrix
     Eigen::Matrix3d  L;
     fstream file;
     std::string word;
     file.open(filename.c_str());
     int count=0;
     int i=0;
     while(file >> word) {
           int j=count%3;
    	   L(i,j)=std::stod(word);
	       count++;
           if(count%3==0){i++;}
           
     }
     return L;
}


void print_ptgroup(const std::vector<Eigen::Matrix3d> point_group, Eigen::Matrix3d lattice)
{
//Prints the number of symmops and then lists each of them
      int numops=point_group.size();
      std::cout<< "Point Group has " << numops << "Symmetry Operations"<<endl;
      Eigen::Vector3d no_translation;
      no_translation << 0.0, 0.0, 0.0;

      int count=1;
      for(Eigen::Matrix3d op:point_group){
          ::Symmetry_Operation symop(op, no_translation);
          for (int i=0; i<3; i++){
              for (int j=0; j<3; j++){
                  symop.cart_matrix(i,j)=op(i,j);
              }
          }
          std::string op_type = check_op_type(symop, lattice);
          std::cout<<"Symmetry Op "<<count<<" is a "<<op_type<<endl;
          std::cout<<op<<endl;
          count++;
      }
      return;
}

void print_factor_group(const std::vector<Symmetry_Operation> factor_group)
{
    int numops=factor_group.size();
    std::cout<< "The Factor Group has " << numops << "Symmetry Operations"<<endl;

    int count=1;
    for(Symmetry_Operation op:factor_group){
          std::cout<<"Symmetry Op " <<count<<endl;
          std::cout<<"Cartesian Matrix"<<endl;
          std::cout<<op.cart_matrix<<endl;
          std::cout<<"Translation Vector"<<endl;
          std::cout<<op.translation<<endl;
          
          count++;
      }
    return;
}

Crystal_Structure read_poscar(std::string filename)
{
    /*input is the path to a Poscar file
     * This function reads the file and puts the information of the lattice and basis
     * into a Crystal_structure object which is returned
     */

     Eigen::Matrix3d  lattice;
     std::vector<Eigen::Vector3d> basis;
     Crystal_Structure xtal_struct(lattice, basis);
     std::ifstream file(filename);
     std::string line;
//     file.open(filename.c_str());
     std::string title;
//     std::getline(file, title);
     int ct=0;
     double val1;
     double val2;
     double val3;
     while (std::getline(file, line)){
             std::stringstream linestream(line);
 //            std::cout<<ct<<std::endl;
             if (ct==0){ct++;
                 title=line;}
             else if (ct==1){ct++; continue;}
             else if (ct>1 && ct<=4){
             linestream>>val1>>val2>>val3;
             xtal_struct.lattice(0,ct-2)=val1;
             xtal_struct.lattice(1,ct-2)=val2;
             xtal_struct.lattice(2,ct-2)=val3;
             ct++;
             }
             else if(ct>4 && ct<8){ct++;}
             else{
             linestream>>val1>>val2>>val3;
             Eigen::Vector3d temp;
             temp(0)=val1;
             temp(1)=val2;
             temp(2)=val3;
             xtal_struct.basis.push_back(temp);
             ct++;
             }
     }
//     std::cout<<"Input file title: "<<title<<std::endl;
//     std::cout<<"Lattice: \n"<<xtal_struct.lattice <<std::endl;
//     std::cout<<"Basis: \n" <<std::endl;
//     int len_basis=xtal_struct.basis.size();
//     for (int i=0; i<len_basis; i++){
//         std::cout<<xtal_struct.basis[i]<<std::endl;
//     }
     return xtal_struct;

}
