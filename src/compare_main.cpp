#include <igl/arap.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>
#include "arap_solver.h"
#include "ply_parser.h"
#include <chrono>

using Clock = std::chrono::steady_clock;
using Milliseconds = std::chrono::duration<double, std::milli>;


int main() {
    //mode = 0 using own arap; mode = 1 using libigl as reference
    int mode = 0;
    
    //Displacement for both arap for comparison
    //Eigen::Vector3d displacement(0.0, 50.0, 0.0);
    //maximal iteration
    int iteration = 20;
    //For visualization
    Eigen::MatrixXd V_display;
    int highlight_handle_index;
    Eigen::MatrixXi F_display;

    if (mode == 0) {
        arap::Mesh mesh = arap::PLYParser::load("data/Simplified_Armadillo.ply");
        Eigen::MatrixXd V0 = mesh.vertices();
        //bounding box diagonal
        Eigen::RowVector3d bbox_min = V0.colwise().minCoeff();
        Eigen::RowVector3d bbox_max = V0.colwise().maxCoeff();
        double diagonal = (bbox_max - bbox_min).norm();
        const Eigen::Vector3d displacement(0.0, 0.1 * diagonal, 0.0);
        //A simple test
        int fixed_index = 0;
        V0.col(1).minCoeff(&fixed_index);
        int handle_index = 0;
        V0.col(1).maxCoeff(&handle_index);

        const int handle = handle_index;
        std::vector<int> fixed = {fixed_index};
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
        const auto time_start = Clock::now();
        solver.initialize(mesh);
        solver.set_constraints(constraint_indices, constraint_positions);
        solver.solve(iteration);
        const auto time_end = Clock::now();

        //calculate time
        double arap_time_ms = Milliseconds(time_end - time_start).count();
        double average_iteration_ms = arap_time_ms / iteration;
        //output
        std::cout << "Solve time for " << iteration << " iterations: " << arap_time_ms << "ms." << std::endl;
        std::cout << "Average iteration time: " << average_iteration_ms << "ms." << std::endl;

        V_display = solver.deformed_vertices();
        highlight_handle_index = handle;
        F_display = mesh.faces();
    }else {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        if(!igl::read_triangle_mesh("data/Simplified_Armadillo.ply", V, F)){
            std::cerr << "Failed to read mesh.\n";
            return 1;
        }
        Eigen::RowVector3d bbox_min = V.colwise().minCoeff();
        Eigen::RowVector3d bbox_max = V.colwise().maxCoeff();
        double diagonal = (bbox_max - bbox_min).norm();
        const Eigen::Vector3d displacement(0.0, 0.1 * diagonal, 0.0);

        int fixed_index = 0;
        V.col(1).minCoeff(&fixed_index);
        int handle_index = 0;
        V.col(1).maxCoeff(&handle_index);

        const int handle = handle_index;
        std::vector<int> fixed_vertices = {fixed_index};
        Eigen::Vector3d target = V.row(handle).transpose() + displacement;
        
        int num_constraints = fixed_vertices.size() + 1;
        Eigen::VectorXi b(num_constraints);
        Eigen::MatrixXd bc(num_constraints, 3);
        for (int i = 0; i < fixed_vertices.size(); ++i) {
            const int vertex_index = fixed_vertices[i];
            b(i) = vertex_index;
            bc.row(i) = V.row(vertex_index);
        }
        b(fixed_vertices.size()) = handle;
        bc.row(fixed_vertices.size()) = target;
        
        igl::ARAPData arap_data;
        arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
        arap_data.max_iter = iteration;
        arap_data.with_dynamics = false;
        
        const auto time_start = Clock::now();
        if (!igl::arap_precomputation(V, F, 3, b, arap_data)) {
            std::cerr << "precomputation failed\n";
            return 1;
        }
        if (!igl::arap_solve(bc, arap_data, V)) {
            std::cerr << "ARAP solve failed\n";
            return 1;
        }
        const auto time_end = Clock::now();

        //calculate time
        double arap_time_ms = Milliseconds(time_end - time_start).count();
        double average_iteration_ms = arap_time_ms / iteration;
        //output
        std::cout << "Solve time(libigl) for " << iteration << " iterations: " << arap_time_ms << "ms." << std::endl;
        std::cout << "Average iteration time(libigl): " << average_iteration_ms << "ms." << std::endl;


        V_display = V;
        highlight_handle_index = handle;
        F_display = F;
    }

    //visualization
    igl::opengl::glfw::Viewer viewer;

    viewer.data().set_mesh(V_display, F_display);
    viewer.data().set_colors(Eigen::RowVector3d(0.7, 0.75, 0.85));
    viewer.data().show_lines = true;

    Eigen::MatrixXd handle_point(1, 3);
    handle_point.row(0) = V_display.row(highlight_handle_index);

    viewer.data().add_points(handle_point, Eigen::RowVector3d(1.0, 0.0, 0.0));
    viewer.data().point_size = 15.0f;

    //set handle point as a center of the camera
    viewer.core().align_camera_center(V_display, F_display);
    auto& core = viewer.core();

    core.orthographic = true;
    core.camera_zoom = 1.0f;
    core.camera_translation = Eigen::Vector3f::Zero();

    constexpr float pi = 3.14159265358979323846f;
    core.trackball_angle =Eigen::Quaternionf(Eigen::AngleAxisf(pi, Eigen::Vector3f::UnitY()));

    viewer.launch();
    return 0;
}