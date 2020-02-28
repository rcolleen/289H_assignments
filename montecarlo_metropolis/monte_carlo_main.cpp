#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "mc_metropolis.cpp"
#include <fstream>
//#include "../matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;
std::vector<double> run_mc_metrolopis(const double& T, const double& mu);

int main()
{
    int Tnum=29, Munum=21;
    double Tmax=1500.0, Tmin=100.0, Mumax=0.50, Mumin=-0.50;
    Eigen::ArrayXd Temps=Eigen::ArrayXd::LinSpaced(Tnum, Tmin, Tmax);
    Eigen::ArrayXd Mus=Eigen::ArrayXd::LinSpaced(Munum, Mumin, Mumax);

    Eigen::MatrixXd Average_energy(Tnum, Munum), Average_xx(Tnum, Munum), Average_Cp(Tnum, Munum);
    Eigen::VectorXd xx =Average_xx.row(1);
    Eigen::MatrixXd configuration = -1.0*Eigen::MatrixXd::Ones(grid_size, grid_size);

    std::ofstream compfile, cpfile;
    compfile.open("MC_composition.txt", std::ofstream::app); cpfile.open("MC_cv.txt", std::ofstream::app);
    compfile<< "Composition wrt Temp (rows) from "<<Tmin<<"to "<<Tmax<<" by Mu (col) from "<<Mumin<<" to "<<Mumax<<".\n";
    cpfile<< "Heat Capacity wrt Temp (rows) from "<<Tmin<<"to "<<Tmax<<" by Mu (col) from "<<Mumin<<" to "<<Mumax<<".\n";
   
    for(int i = 0; i< Munum; i++)
    { 
        double Mu =Mus(i);
        for (int j=0; j< Tnum; j++)
        {
            double T = Temps(j);
            std::vector<double> temp = run_mc_metropolis(T, Mu, &configuration);
            std::cout<<"Metropolis Output at Mu =" <<Mu<<": " <<std::endl;
            std::cout<<"Metropolis Output at Temp =" <<T<<": " <<std::endl;
            std::cout<< temp[0]<< "\t"<<temp[1]<<std::endl;
//          Average_energy(i, j) = temp[0];
            Average_xx(j,i) = temp[0];
            Average_Cp(j,i) = temp[1];
            cpfile<<Average_Cp(j,i)<<"\t";

        }
        cpfile<<"\n";
        std::string ttl ="Mu: " +std::to_string(Mu);
        plt::title(ttl);
        plt::xlabel("Atomic Fraction A");
        plt::ylabel("Temperature (K)");
//        std::cout<<"Average_xx at T = "<<T<<": \n"<<Average_xx.row(i)<<std::endl;
//        std::cin.get();
        xx =Average_xx.col(i).transpose();
        for (int k = 0; k<Tnum;k++){
            xx(k)=xx(k)/(grid_size*grid_size);
            Average_xx(k,i)=xx(k);
            compfile<<Average_xx(k,i)<<"\t";
        }
//        std::cout<<"Print XX: \n"<<Average_xx<<std::endl;
        compfile<<"\n";

        plt::plot(xx, Temps, "*r");
//      std::cin.get();

//      plt::show();
        plt::pause(0.01);

    }
    cpfile.close(); compfile.close();
 
    plt::figure(2);
    for(int k=0; k<Tnum;k++)
    {
        plt::plot(Average_xx.row(k), Mus, "*");
    }
    plt::xlabel("Atomic Fraction A");
    plt::ylabel("Chemical Potential");
    plt::show();
//    std::cin.get();

    return 0;
}


