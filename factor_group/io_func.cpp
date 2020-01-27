#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <set>
//#include "symop.cpp"
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

const char *symop_label_table[48] = {
    "A",
    "B",
    "C",
    "D",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
    "AA",
    "AB",
    "AC",
    "AD",
    "AE",
    "AF",
    "AG",
    "AH",
    "AI",
    "AJ",
    "AK",
    "AL",
    "AM",
    "AN",
    "AO",
    "AP",
    "AQ",
    "AR",
    "AS",
    "AT",
    "AU",
    "AV",
    "AW"
};



//auto find_multiplication_table(const std::vector<Eigen::Matrix3d>& point_group)
//{
//    int grp_size = point_group.size();
//    char multiplication_table[grp_size][grp_size];
//
//    return multiplication_table;
//}
//
//
//void print_multiplication_table(const std::vector<Eigen::Matrix3d>& point_group)
//{
////determines the multiplcation table for the provided point group and 
//// then prints it to the screen
//    int grp_size = point_group.size();
//    char multiplication_table[grp_size][grp_size] = find_multiplication_table(point_group);
//    std::cout<< "Multiplication Table for Point Group:" <<std::endl;
//    std::cout<<multiplication_table<<std::endl;
//    return;
//}
//
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

