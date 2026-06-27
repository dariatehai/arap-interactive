#include "ply_parser.h"
#include "arap_solver.h"
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

#include <iostream>
#include <stdexcept>

int main()
{
    /*
    try {
        arap::Mesh mesh = arap::PLYParser::load("data/airplane.ply");

        if (!mesh.is_valid()) {
            std::cerr << "Error: Loaded mesh is empty." << std::endl;
            return 1;
        }

        const Eigen::MatrixXd& V = mesh.vertices();
        const Eigen::MatrixXi& F = mesh.faces();

        std::cout << "Loaded mesh: "
                  << V.rows() << " vertices, "
                  << F.rows() << " faces." << std::endl;

        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(V, F);
        viewer.data().set_face_based(true);
        viewer.launch();

    } catch (const std::exception& e) {
        std::cerr << "Error loading PLY file: " << e.what() << std::endl;
        return 1;
    }
    */

    try{
        //--------------------------Load mesh-----------------------------------
        arap::Mesh mesh = arap::PLYParser::load("data/Simplified_Armadillo.ply");
        if (!mesh.is_valid()) {
            std::cerr << "Error: Loaded mesh is empty." << std::endl;
            return 1;
        }
        const Eigen::MatrixXd& V0 = mesh.vertices();
        const Eigen::MatrixXi& F = mesh.faces();
        
        //-------------------Initialize arapSolver--------------------------------
        arap::ARAPSolver solver;
        solver.initialize(mesh);

        //----------------simple constraints(only 2 constraints)---------------------------
        int fixed_vertex = 0;
        int handle_vertex = static_cast<int>(V0.rows()) / 2;

        std::vector<int> constraint_indices = {
            fixed_vertex,
            handle_vertex
        };

        Eigen::MatrixXd constraint_positions(2, 3);

        // fixed vertex stays unchanged
        constraint_positions.row(0) = V0.row(fixed_vertex);

        // handle vertex is moved
        constraint_positions.row(1) = V0.row(handle_vertex);

        Eigen::RowVector3d min_v = V0.colwise().minCoeff();
        Eigen::RowVector3d max_v = V0.colwise().maxCoeff();
        Eigen::RowVector3d extent = max_v - min_v;

        double move_scale = 0.25 * extent.norm();

        constraint_positions.row(1) += Eigen::RowVector3d(
            move_scale,
            0.0,
            0.0
        );

        //--------------------------------solve arap-----------------------------------------
        solver.set_constraints(constraint_indices, constraint_positions);
        solver.solve(10);

        Eigen::MatrixXd V_deformed = solver.deformed_vertices();
        //----------------------------show mesh side by side------------------------------------------------
        int n = static_cast<int>(V0.rows());
        int m = static_cast<int>(F.rows());

        Eigen::MatrixXd V_original_show = V0;
        Eigen::MatrixXd V_deformed_show = V_deformed;

        double shift = 1.0 * extent.norm();

        V_original_show.col(0).array() -= shift;
        V_deformed_show.col(0).array() += shift;

        Eigen::MatrixXd V_show(2 * n, 3);
        V_show.topRows(n) = V_original_show;
        V_show.bottomRows(n) = V_deformed_show;

        Eigen::MatrixXi F_show(2 * m, 3);
        F_show.topRows(m) = F;
        F_show.bottomRows(m) = F.array() + n;

        Eigen::MatrixXd colors(2 * m, 3);

        // original mesh: gray
        colors.topRows(m).rowwise() = Eigen::RowVector3d(0.75, 0.75, 0.75);

        // deformed mesh: orange
        colors.bottomRows(m).rowwise() = Eigen::RowVector3d(1.0, 0.55, 0.15);

        igl::opengl::glfw::Viewer viewer;

        viewer.data().set_mesh(V_show, F_show);
        viewer.data().set_colors(colors);
        viewer.data().show_lines = true;
        //-------------------------add constaint points visually--------------------------------
        Eigen::MatrixXd points(2, 3);

        points.row(0) = constraint_positions.row(0);
        points.row(1) = constraint_positions.row(1);

        // shift points to the deformed mesh side
        points.col(0).array() += shift;

        Eigen::MatrixXd point_colors(2, 3);

        // blue = fixed vertex
        point_colors.row(0) = Eigen::RowVector3d(0.0, 0.0, 1.0);

        // red = moved handle vertex
        point_colors.row(1) = Eigen::RowVector3d(1.0, 0.0, 0.0);

        viewer.data().add_points(points, point_colors);

        viewer.launch();

    }catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}