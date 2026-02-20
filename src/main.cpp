// src/main.cpp
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/TriangularMesh.h"
#include </mnt/c/Users/mmaru/Desktop/AE623/Project 2/external/eigen/Eigen/Dense>
#include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/bj_limiter.h"
#include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/StateMesh.h"


int main() {
    std::shared_ptr<TriangularMesh> mesh = std::make_shared<TriangularMesh>("projects/Project-1/mesh_refined_2394.gri");

    // Inlet conditions
    double gamma = 1.4;
    double rho0 = 1;
    double a0 = 1;
    double p0 = rho0*a0*a0/gamma;
    double alpha = 50*3.14/180;
    double pout = 0.7*p0;
    double M = 0.1;

    // Boundary conditions
    std::vector<std::shared_ptr<BoundaryCondition>> bc{};

    // Initialize the state mesh
    StateMesh states(mesh, bc);

    // auto mesh = std::make_shared<TriangularMesh>("projects/Project-1/mesh_coarse.gri");
    // // TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");
    // // auto meshPtr = std::make_shared<mesh>;
    // StateMesh states(mesh);
    std::vector<Eigen::Matrix<double,4,2>> Lgrad(states.cellCount(), Eigen::Matrix<double,4,2>::Constant(0.5));
    // Eigen::Matrix3d Lgrad =  Eigen::Matrix3d::Constant(0.5);
    std::vector<Eigen::Matrix<double,4,2>> Lgrad_new;
    Lgrad_new = bj_limiter(Lgrad,states);

    // mesh.writeGri("projects/Project-1/mesh_coarse.gri");

    return 0;
}