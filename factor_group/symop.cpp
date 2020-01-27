#include "Eigen/Core"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>

class Symmetry_Operation
{

public:
    Eigen::Matrix3d cart_matrix;
    Eigen::Vector3d translation;

    Symmetry_Operation(Eigen::Matrix3d input_matrix,  Eigen::Vector3d input_translation)
    {
        Eigen::Matrix3d cart_matrix;
        for(int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                cart_matrix(i,j)=input_matrix(i,j);
            }}
//        std::cout<<"CART Matrix: "<< cart_matrix<<std::endl;
        Eigen::Vector3d translation = input_translation;
    }
};
