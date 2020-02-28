#include "Eigen/Dense"
#include <random>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <math.h>
#include "../matplotlib-cpp/matplotlibcpp.h"

#define PREC 1e-6
namespace plt = matplotlibcpp;

//Global Variables declared (for monte carlo parameters)
//
const int grid_size=100;
const double Vnn=0.030, kb=8.617333e-5, Vpt=0.0, V0=-2*Vnn;

int periodic_assign(int index, int length)
{
    return (length+(index%length))%length;
}

double random_number()
{
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    return dis(gen);
}

double calc_delta_E(const Eigen::MatrixXd& configuration, int index_i, int  index_j)
{
    int index_ip=periodic_assign(index_i+1, grid_size), index_im=periodic_assign(index_i-1 ,grid_size); 
    int index_jp=periodic_assign(index_j+1, grid_size), index_jm=periodic_assign(index_j-1, grid_size); 
//    std::cout<<"Neighbors: i+1="<<index_ip<<", i-1= "<<index_im<<", j+1= "<<index_jp<<", j-1= "<<index_jm<<std::endl;
    double delta_E = -2*(Vpt*configuration(index_i, index_j)+Vnn*configuration(index_i,index_j)*(configuration(index_i, index_jp)+configuration(index_i, index_jm)+configuration(index_ip,index_j)+configuration(index_im,index_j)));
    return delta_E;
}

double calc_delta_grand_canonical_energy(const Eigen::MatrixXd& configuration, const int& index_i, const int& index_j, const double& chem_potential)
{
    double delta_E = calc_delta_E(configuration, index_i, index_j);
    return delta_E-((-1.0*configuration(index_i,index_j))*chem_potential);
}

double calc_total_N(Eigen::MatrixXd configuration)
{
    Eigen::MatrixXd temp= configuration;
    temp = configuration.array()+1.0;
    return temp.sum()/2;
}

double calc_total_grand_canonical_energy(const Eigen::MatrixXd& configuration, const double& chemical_potential)
{
    //solving (E(sigma)-N*mu)
    double chem_potential_total;
    double N= calc_total_N(configuration);
//    std::cout<<"Total N = "<<N<<std::endl;
    chem_potential_total = N*chemical_potential;
    double config_energy= V0+configuration.sum()*Vpt;
  //  std::cout<<"Configuration.sum() = "<<configuration.sum()<<std::endl;
    for (int i=0; i<grid_size; i++){
        for(int j=0; j<grid_size; j++){
            double temp=Vnn*configuration(i,j)*(configuration(i,((j+1)%grid_size))+ configuration(((i+1)%grid_size),j));
            config_energy+=temp;
        }
    }
    double energy;
    energy= (config_energy-chem_potential_total);
    return energy;
}



bool accept_change(const Eigen::MatrixXd& configuration, const int& index_i, const int& index_j, const double& temperature, const double& chem_potential, std::mt19937_64& generator, std::uniform_real_distribution<>& distribution)//, std::vector<int> *rate_vect_ptr)
{
    /* This function calculated the delta energy. If negative it accepts. 
     * If exp(energy/kbT)> rand number then it is also accepted
     * Otherwise reject.
     */
    double delta_gce = calc_delta_grand_canonical_energy(configuration,index_i, index_j, chem_potential);
    double en_test = std::exp(-1.0*(delta_gce)/(kb*temperature));
//    std::vector<int>& rate_vect= *rate_vect_ptr;
    if (delta_gce<0){
 //       std::cout<< "Accept Condition 1: delta_gce = "<<delta_gce<<std::endl;
//        rate_vect.push_back(0);
        return true;
    }
    else if (en_test>=distribution(generator)){
 //       std::cout<< "Accept Condition 2: delta_gce = "<< delta_gce<<std::endl;
  //      rate_vect.push_back(1);
        return true;
    }
    else{
 //       std::cout<<"Not Accepted"<<std::endl; 
    ///    rate_vect.push_back(2);
        return false;
    }
}


double find_vector_average(const std::vector<double>& vect)
{
    double vsum = 0.0;
//    double ct=0;
    for(double val:vect){
        vsum+=val;
  //      ct++;
    }
    double vsize=vect.size();
    vsum = vsum/vsize;
 //   std::cout<<"vector average in function: "<<vsum<<std::endl;
    return vsum;

}

double find_cp_average(const std::vector<double>& energy,double T)
{
    double e_average_sq=find_vector_average(energy);
    std::cout<<"average_energy = "<<e_average_sq<<std::endl;
    e_average_sq = pow(e_average_sq,2.0);
    auto e_squares = energy;

    for(int i=0; i<int(energy.size()); i++)
    {   //std::cout<<"e = "<<energy[i];
        e_squares[i]=pow(e_squares[i],2.0);
        //std::cout<<"; e^2 = "<<e_squares[i]<<std::endl;
    }
    double average_of_esquares = find_vector_average(e_squares);
    double kbt2=kb*pow(T,2.0);
    double cv = average_of_esquares-e_average_sq;
    cv=cv/kbt2;
    std::cout<<"number of energies = "<<energy.size()<<std::endl;
    std::cout<<"average of squared energies = "<<average_of_esquares<<std::endl;
    std::cout<<"square of average energy = "<<e_average_sq<<std::endl;
    std::cout<<"kb*T^2 = "<<kbt2<<std::endl;
    return cv;

}

void calc_accept_rates(std::vector<int> rates_vect)
{
    int accept1=0, accept2=0, reject=0;
    for(int val:rates_vect)
    {
        if (val==0){accept1++;}
        else if (val==1){accept2++;}
        else if (val==2){reject++;}
    }
    double accept1_rate=double(accept1)/rates_vect.size(), accept2_rate=double(accept2)/rates_vect.size(), reject_rate=double(reject)/rates_vect.size();

    std::cout<<"Acceptance rate:\n";
    std::cout<<"Accept1 = "<<accept1_rate<<std::endl;
    std::cout<<"Accept2 = "<<accept2_rate<<std::endl;
    std::cout<<"Reject  = "<<reject_rate<<std::endl;
}

std::vector<double> run_mc_metropolis(const double& temperature, const double& chem_potential, Eigen::MatrixXd *configuration_pointer)
{
// generates a nxN grid with a Vpt and Vnn ECI's set and uses a MC metropolis algorithm
// to find the average grand canonical energy for a given T and chemical potential

    Eigen::MatrixXd& configuration = *configuration_pointer;
    int max_pass = 7000, index_i, index_j, num_per_pass=grid_size*grid_size, thresh=2000;
//    double  delta_average_energy=1.0;//threshold = 1e-4,
    std::vector<double> average_energy, average_N, averages;
 //   average_energy.reserve(max_pass-thresh);
 //   average_N.reserve(max_pass-thresh);

//    std::cout<<"Initital Configuration"<< configuration<<std::endl;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1

    double total_energy=calc_total_grand_canonical_energy(configuration, chem_potential);
    double total_N =calc_total_N(configuration);
//    std::vector<int> accept_rate;
    for(int pass_num=0;pass_num<max_pass+thresh; pass_num++)
    {
         for(int loop_in_pass=0; loop_in_pass<num_per_pass; loop_in_pass++)
         {
             index_i=int(dis(gen)*grid_size) % grid_size;
             index_j=int(dis(gen)*grid_size) % grid_size; 
             double L_ij=-1*configuration(index_i, index_j);  
             if(accept_change(configuration, index_i, index_j, temperature, chem_potential, gen, dis))//accept_rate))
             {
                 total_energy+=calc_delta_grand_canonical_energy(configuration,index_i, index_j, chem_potential);
                 total_N+=L_ij;
                 configuration(index_i, index_j)=L_ij;
             }
         }
         if (pass_num>=thresh)
         {
             average_N.push_back(total_N);
             average_energy.push_back(total_energy);
//             std::cout<<"total N = "<< total_N<<std::endl;
//            std::cout<<"total_energy = "<<total_energy<<std::endl;
         }
    }
    averages.push_back(find_vector_average(average_N));
    averages.push_back(find_cp_average(average_energy, temperature));
//    std::cout<<"Configuration: \n"<<configuration<<std::endl;
//    calc_accept_rates(accept_rate);
    return averages;
}
