#include "Eigen/Dense"
#include <iostream>
#include <cmath>
#include <fstream>
#include "io_func.cpp"
#include "point_group.cpp"
#include <vector>

#define PREC 1e-6

using namespace std;

int main(int argc, char *argv[])
{ 
//main function which reads in lattice from input file and calcualted point group
//and prints it to the screen
    if (argc != 2){
        std::cout<<"WRONG NUMBER OF INPUT ARGUMENTS!"<<std::endl;
        std::cout<<"Please specify path to POSCAR file after executable."<<std::endl;
        return 1;
    }
    Eigen::Matrix3d L = read_lattice(argv[1]);
    std::cout<<"This is the Lattice:"<<endl;
    std::cout<<L<<endl;
    std::vector<Eigen::Matrix3d> pt_grp;
    bool closed;

    pt_grp = calc_point_group(L);
    closed = is_closed(pt_grp);
    std::cout<<"This Group is Closed: "<<closed<<endl;
//    std::cout<<pt_grp.size()<<endl;
    print_ptgroup(pt_grp);
    return 0;
}
