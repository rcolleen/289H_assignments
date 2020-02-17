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
int grid_size=8;
float Vnn=0.030, kb=8.617333e-5, Vpt=0.0, V0=-2*Vnn;

int periodic_assign(int index, int length)
{
    return (length+(index%length))%length;
}

float random_number()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
    return dis(gen);
}

float calc_delta_E(Eigen::MatrixXf configuration, int index_i, int  index_j)
{
    int index_ip=periodic_assign(index_i+1, grid_size), index_im=periodic_assign(index_i-1 ,grid_size); 
    int index_jp=periodic_assign(index_j+1, grid_size), index_jm=periodic_assign(index_j-1, grid_size); 
//    std::cout<<"Neighbors: i+1="<<index_ip<<", i-1= "<<index_im<<", j+1= "<<index_jp<<", j-1= "<<index_jm<<std::endl;
    float delta_E = -2*(Vpt*configuration(index_i, index_j)+Vnn*configuration(index_i,index_j)*(configuration(index_i, index_jp)+configuration(index_i, index_jm)+configuration(index_ip,index_j)+configuration(index_im,index_j)));
    return delta_E;
}

float calc_delta_grand_canonical_energy(Eigen::MatrixXf configuration, int index_i, int index_j, const float& chem_potential)
{
    float delta_E = calc_delta_E(configuration, index_i, index_j);
    return delta_E-((-1.0*configuration(index_i,index_j))*chem_potential);
}

float calc_total_N(Eigen::MatrixXf configuration)
{
    Eigen::MatrixXf temp = configuration;
    temp=temp.array()+1.0;
//    for (int i=0; i<temp.size(); i++)
  //      { temp(i)++;}
    return temp.sum()/2;
}

float calc_total_grand_canonical_energy(Eigen::MatrixXf configuration, const float& chemical_potential, const float& temperature)
{
    //solving (E(sigma)-N*mu)/kbT
    float chem_potential_total;
    chem_potential_total = calc_total_N(configuration)*chemical_potential;
    float config_energy= V0+configuration.sum()*Vpt;
    for (int i=0; i<grid_size; i++){
        for(int j=0; j<grid_size; j++){
            float temp=Vnn*configuration(i,j)*(configuration(i,((j+1)%grid_size))+ configuration(((i+1)%grid_size),j));
            config_energy+=temp;
        }
    }
    float energy;
    energy= (config_energy-chem_potential_total)/(kb*temperature);
    return energy;
}



bool accept_change(Eigen::MatrixXf configuration, int index_i, int index_j, const float& temperature, const float& chem_potential)
{
    /* This function calculated the delta energy. If negative it accepts. 
     * If exp(energy/kbT)> rand number then it is also accepted
     * Otherwise reject.
     */
    float delta_gce = calc_delta_grand_canonical_energy(configuration,index_i, index_j, chem_potential);
    float en_test = std::exp(-1.0*(delta_gce)/(kb*temperature));
    if (delta_gce<0){
 //       std::cout<< "Accept Condition 1: delta_gce = "<<delta_gce<<std::endl;
        return true;
    }
    else if (en_test>random_number()){
 //       std::cout<< "Accept Condition 2: delta_gce = "<< delta_gce<<std::endl;
        return true;
    }
    else{
 //       std::cout<<"Not Accepted"<<std::endl; 
        return false;
    }
}

float update_average_quantity( float& old_average, float& new_quantity, int& pass_num)
{
    /* This function takes a running average a new quantity and number of the current pass and 
     * uses this info to update the running average.
     */
    float new_average;
    if (pass_num>1){
        new_average = old_average*float(pass_num/(pass_num+1.0))+new_quantity*float(1.0/(pass_num+1.0));}
    else{ new_average = new_quantity;}
    return new_average;

}

float find_vector_average(std::vector<float>& vect)
{
    float vsum = 0.0;
//    float ct=0;
    for(float val:vect){
        vsum+=val;
  //      ct++;
    }
    float vsize=vect.size();
    vsum = vsum/vsize;
 //   std::cout<<"vector average in function: "<<vsum<<std::endl;
    return vsum;

}

float find_cp_average(std::vector<float>& energy)
{
    float esum=0.0;
    float esq_sum=0.0;
    float energy_sz=energy.size();
    for(int i=0; i<int(energy.size()); i++)
    {
        esum+=energy[i];
        esq_sum+=pow(energy[i],2.0);
    }
    return (pow(esum/energy_sz,2.0)+esq_sum)/energy_sz;

}



std::vector<float> run_mc_metropolis(const float& temperature, const float& chem_potential)
{
// generates a nxN grid with a Vpt and Vnn ECI's set and uses a MC metropolis algorithm
// to find the average grand canonical energy for a given T and chemical potential

    Eigen::MatrixXf configuration = -1.0*Eigen::MatrixXf::Ones(grid_size,grid_size);
    int min_i = 50, max_i = 5000, pass_num=1, index_i, index_j;
    float  delta_average_energy=1.0;//threshold = 1e-4,
    std::vector<float> average_energy, average_N, averages;

    float new_energy= calc_total_grand_canonical_energy(configuration, chem_potential, temperature);
    float new_average_N= calc_total_N(configuration);
    average_N.push_back(new_average_N);
//    std::cout<<"Initital Configuration"<< configuration<<std::endl;
    while (pass_num<max_i){
        index_i=int(random_number()*grid_size);
        index_j=int(random_number()*grid_size);
  //      std::cout<<"Index I and J: " <<index_i<<", "<<index_j<<std::endl;
        float L_ij=configuration(index_i, index_j);
   //     std::cout<<"L_ij = "<<L_ij<<std::endl;
        
        if(accept_change(configuration, index_i, index_j, temperature, chem_potential)){
           configuration(index_i, index_j)=-1*L_ij;}
     
     //      std::cout <<"Accepted config# "<<pass_num<<std::endl;
        new_energy= calc_total_grand_canonical_energy(configuration, chem_potential, temperature);
        new_average_N= calc_total_N(configuration);
//           std::cout<<"Total N = "<< new_average_N<<std::endl;
        average_N.push_back(new_average_N);
    //       std::cout<<"new energy :"<<new_energy<<std::endl;
    //       if (pass_num>min_i){
    //         float prev_average_e=find_vector_average(average_energy);
    //           std::cout<<"Average_energy:  "<<prev_average_e<<std::endl;
        average_energy.push_back(new_energy);
        //       float new_average_e=find_vector_average(average_energy);
         //      delta_average_energy=std::abs(prev_average_e-new_average_e);
   //            std::cout<<"Delta Average Energy: "<<delta_average_energy;
          //   average_energy.push_back(new_energy);
      //     }
     //      else average_energy.push_back(new_energy);
        pass_num++; //only increments when a new config is accepted
        }
   //     if (pass_num>max_i){
   //            break;
   //        }}
    
    
//    std::cout<<"Final Pass Number:  "<<pass_num<<std::endl;
//    std::cout<<"Delta Average Energy: "<<delta_average_energy;
//    std::vector<float> averages;
//    averages.push_back(find_vector_average(average_energy));
//    float average_xx = find_vector_average(average_N)/float(grid_size^2);
    averages.push_back(find_vector_average(average_N));
//    std::cout<<"Average_N: "<<find_vector_average(average_N)<<std::endl;
//    plt::hist(average_N, 10);
//    plt::show();
//    std::cin.get();
    averages.push_back(find_cp_average(average_energy));
    return averages;
}
