#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "Eigen/src/Core/Matrix.h"
#include "Eigen/src/Core/util/Constants.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "symop.cpp"
#include "crystal_class.cpp"
#include "point_group.cpp"

#define PREC 1e-6


auto fractional_to_cartesian(const Eigen::Matrix3d lattice, Eigen::Vector3d coordinate)
{
    Eigen::Vector3d cart_coordinate;
    cart_coordinate(0)=coordinate(0)*lattice(0,0)+coordinate(1)*lattice(0,1)+coordinate(2)*lattice(0,2);
    cart_coordinate(1)=coordinate(1)*lattice(1,0)+coordinate(1)*lattice(1,1)+coordinate(2)*lattice(1,2);
    cart_coordinate(2)=coordinate(2)*lattice(2,0)+coordinate(1)*lattice(2,1)+coordinate(2)*lattice(2,2);
    return cart_coordinate;
}

auto apply_symop_to_basis(Symmetry_Operation this_op, std::vector<Eigen::Vector3d> basis)
{
    /* this steps through each coordinate in the basis and applies the given symmetry 
     * operation to it. The final vector of transformed coordinates is returned.
     */
    std::vector<Eigen::Vector3d> basis_prime;
    for(const Eigen::Vector3d& coord : basis){
        Eigen::Vector3d coord_prime;
        coord_prime=this_op.cart_matrix*coord+this_op.translation;
        basis_prime.push_back(coord_prime);
    }
    return basis_prime;
}

struct vector_is_same{
    vector_is_same( Eigen::Vector3d vect_a) : vect_a(vect_a) {}
    bool operator()( Eigen::Vector3d vect_b){
        if(vect_a.isApprox(vect_b,PREC) ){return true;}
        else return false;
    }
private:
    Eigen::Vector3d vect_a;
};

bool maps_to_original_basis(std::vector<Eigen::Vector3d> new_basis, std::vector<Eigen::Vector3d> original_basis)
{
    /*will iterate through each new basis member and check if it also exists in the original basis
     *if one element is found that does not exist in the original basis, return false
     *if all elements are found, return true
     */
    int len=new_basis.size();
    for (int i=0; i<len;i++){
        Eigen::Vector3d coord_prime=new_basis[i];
        vector_is_same coord_compare(coord_prime);
        if(std::find_if(original_basis.begin(), original_basis.end(), coord_compare)==original_basis.end()){
            return false;
        }
    }
    return true;

}

auto translate_new_basis(std::vector<Eigen::Vector3d> transformed_basis, Eigen::Vector3d translation)
{
    /*applied translation vector to every transformed basis element and returns new basis
     *
     */
    int len=transformed_basis.size();
    std::vector<Eigen::Vector3d> new_basis;
    for (int i=0; i<len; i++){
        Eigen::Vector3d coord_original= transformed_basis[i];
        auto coord_prime=coord_original+translation;
//        coord_prime= map_to_unit_cell(coord_prime)
        new_basis.push_back(coord_prime);
    }
    return new_basis;
}

std::vector<Symmetry_Operation> find_factor_group(Crystal_Structure xtal_struct)
{
/*takes crystal structure as input (lattice and basis), finds point group, then checks each point 
 *group operation, with all translations, to see if it remaps to the original basis. if it does,
 *the symop is added to the factor group. The factor group is returned.
 */
    auto lattice = xtal_struct.lattice;
    auto basis = xtal_struct.basis;

    auto pt_group = calc_point_group(lattice);
    int len_group = pt_group.size();
    int len_basis = basis.size();
    Eigen::Matrix3d dummy_cart = Eigen::Matrix3d::Zero();
    Eigen::Vector3d dummy_trans = Eigen::Vector3d::Zero();
    Symmetry_Operation this_symop(dummy_cart, dummy_trans);
    
    std::vector<Symmetry_Operation> factor_group;

    for(int i=0; i<len_group; i++){
         for(int j = 0; j<3; j++){
             for(int k=0; k<3;k++){
                 this_symop.cart_matrix(j,k)=pt_group[i](j,k);
             }
         }
         this_symop.translation<< 0,0,0;
         auto transformed_basis = apply_symop_to_basis(this_symop, basis);
         if(maps_to_original_basis(transformed_basis, basis)){factor_group.push_back(this_symop);}
         else
            for (int r1=0; r1<len_basis; r1++){
                for(int r2=0; r2<len_basis; r2++){
                    Eigen::Vector3d r_init=basis[r1];
                    Eigen::Vector3d r_prime=transformed_basis[r2];
                    Eigen::Vector3d translation = r_init - r_prime;
                    auto translated_basis = translate_new_basis(transformed_basis, translation);
                    if(maps_to_original_basis(translated_basis,basis)){
                        this_symop.translation = translation;
                        factor_group.push_back(this_symop);
                    }
                }
            }

 
    }
    return factor_group;
}
