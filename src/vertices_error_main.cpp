#include <igl/arap.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include "arap_solver.h"
#include "ply_parser.h"

double compute_relative_rms_error(
    const Eigen::MatrixXd& V0,
    const Eigen::MatrixXd& V_ours,
    const Eigen::MatrixXd& V_libigl)
{
    const Eigen::MatrixXd difference = V_ours - V_libigl;

    const double rms_error = std::sqrt(difference.squaredNorm()/ static_cast<double>(V0.rows()));

    const Eigen::RowVector3d bbox_min = V0.colwise().minCoeff();
    const Eigen::RowVector3d bbox_max = V0.colwise().maxCoeff();
    const double bbox_diagonal = (bbox_max - bbox_min).norm();

    if (bbox_diagonal < 1e-12) {
        throw std::runtime_error("Original mesh has a zero bounding-box diagonal.");
    }

    return rms_error / bbox_diagonal;
}

int main() {
    //Displacement for both arap for comparison
    Eigen::Vector3d displacement(0.0, 50.0, 0.0);
    //maximal iteration
    int iteration = 100;

    arap::Mesh mesh = arap::PLYParser::load("data/Simplified_Armadillo.ply");
    Eigen::MatrixXd V0 = mesh.vertices();
    //A simple test
    const int handle = V0.rows() / 16;
    std::vector<int> fixed = {0};
    Eigen::Vector3d target = V0.row(handle).transpose() + displacement;

    //set indices of constraint
    std::vector<int> constraint_indices = fixed;
    constraint_indices.push_back(handle);
    //set positions of constraint
    Eigen::MatrixXd constraint_positions(constraint_indices.size(), 3);
    //fixed
    for (int row = 0; row < static_cast<int>(fixed.size()); ++row) {
        constraint_positions.row(row) = V0.row(fixed[row]);
    }
    //handle
    constraint_positions.row(constraint_indices.size() - 1) = target;

    arap::ARAPSolver solver;
    solver.initialize(mesh);
    solver.set_constraints(constraint_indices, constraint_positions);
    solver.solve(iteration);

    Eigen::MatrixXd V_ours = solver.deformed_vertices();

    //libigl as reference
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if(!igl::read_triangle_mesh("data/Simplified_Armadillo.ply", V, F)){
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }
        
    int num_constraints = fixed.size() + 1;
    Eigen::VectorXi b(num_constraints);
    Eigen::MatrixXd bc(num_constraints, 3);
    for (int i = 0; i < fixed.size(); ++i) {
        const int vertex_index = fixed[i];
        b(i) = vertex_index;
        bc.row(i) = V.row(vertex_index);
    }
    b(fixed.size()) = handle;
    bc.row(fixed.size()) = target;
        
    igl::ARAPData arap_data;
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
    arap_data.max_iter = iteration;
    arap_data.with_dynamics = false;
        
    
    if (!igl::arap_precomputation(V, F, 3, b, arap_data)) {
        std::cerr << "precomputation failed\n";
        return 1;
    }
    if (!igl::arap_solve(bc, arap_data, V)) {
        std::cerr << "ARAP solve failed\n";
        return 1;
    }

    Eigen::MatrixXd V_libigl = V;
    
    double relative_rms_error = compute_relative_rms_error(V0, V_ours, V_libigl);

    std::cout << "Relative RMS vertex error: " << relative_rms_error * 100.0 << "%\n";
    
    return 0;
}