#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "mc_metropolis.cpp"
//#include "../matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;
std::vector<float> run_mc_metrolopis(const float& T, const float& mu);

int main()
{
    int Tnum=21, Munum=75;
    float Tmax=1000.0, Tmin=100.0, Mumax=0.50, Mumin=-0.50;
    Eigen::ArrayXf Temps=Eigen::ArrayXf::LinSpaced(Tnum, Tmin, Tmax);
    Eigen::ArrayXf Mus=Eigen::ArrayXf::LinSpaced(Munum, Mumin, Mumax);

    Eigen::MatrixXf Average_energy(Tnum, Munum), Average_xx(Tnum, Munum), Average_Cp(Tnum, Munum);
    Eigen::VectorXf xx =Average_xx.row(1);

    for(int i = 0; i< Tnum; i++)
    { 
        float T = Temps(i);
        for (int j=0; j< Munum; j++)
        {
            float Mu =Mus(j);
            std::vector<float> temp = run_mc_metropolis(T, Mu);
          //  std::cout<<"Metrolopis Output at Mu =" <<Mu<<": " <<std::endl;
          //  std::cout<< temp[0]<< "\t"<<temp[1]<<std::endl;
//            Average_energy(i, j) = temp[0];
            Average_xx(i, j) = temp[0];
            Average_Cp(i, j) = temp[1];

        }
        std::string ttl ="Temp: " +std::to_string(T);
        plt::title(ttl);
        plt::xlabel("Atomic Fraction A");
        plt::ylabel("Chemical Potential");
//        std::cout<<"Average_xx at T = "<<T<<": \n"<<Average_xx.row(i)<<std::endl;
//        std::cin.get();
        xx =Average_xx.row(i);
        for (int k = 0; k<Munum;k++){
            xx(k)=xx(k)/(grid_size*grid_size);
            Average_xx(i,k)=xx(k);
        }
//        std::cout<<"Print XX: \n"<<Average_xx<<std::endl;

        plt::plot(xx, Mus, "*");
//      std::cin.get();

//      plt::show();
        plt::pause(0.01);

    }
 
    plt::figure(2);
    for(int k=0; k<Munum;k++)
    {
        plt::plot(Average_xx.col(k).transpose(), Temps, "*r-");
    }
    plt::xlabel("Atomic Fraction A");
    plt::ylabel("Temperature (K)");
    plt::show();
    std::cin.get();

    return 0;
}


