#include "Eigen/Core"
#include "Eigen/src/Core/Matrix.h"
#include <iostream>
#include <vector>

class Crystal_Structure
{
    public:
        Eigen::Matrix3d lattice;
        std::vector<Eigen::Vector3d> basis;

        Crystal_Structure(Eigen::Matrix3d input_matrix, std::vector<Eigen::Vector3d> input_basis)
        { 
            Eigen::Matrix3d lattice = input_matrix;
            std::vector<Eigen::Vector3d> basis = input_basis;
        }
};
