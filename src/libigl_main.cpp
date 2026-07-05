#include <igl/arap.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //read mesh
    if(!igl::read_triangle_mesh("data/Simplified_Armadillo.ply", V, F)){
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }
    //number of constraints
    int num = 2;
    
    //Define constraints b(n): n constraint indices
    Eigen::VectorXi b(num);
    //same indices as in main.cpp
    int handle_vertex = static_cast<int>(V.rows()) / 2;
    b(0) = 0;
    b(1) = handle_vertex;
    //Define target position as a n*3 matrix
    Eigen::MatrixXd bc(num, 3);

    bc.row(0) = V.row(0);   //fixed vertex
    bc.row(1) = V.row(handle_vertex);   //handle vertex

    //To compare with our main, just use the same moving scale as in the main.cpp
    Eigen::RowVector3d min_v = V.colwise().minCoeff();
    Eigen::RowVector3d max_v = V.colwise().maxCoeff();
    Eigen::RowVector3d extent = max_v - min_v;
    double move_scale = 0.25 * extent.norm();
    bc.row(1) += Eigen::RowVector3d(
        move_scale,
        0.0,
        0.0
    );
    //ARAP refering to Github libigl/include/igl/arap.h
    igl::ARAPData arap_data;
    //configure arap
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
    arap_data.max_iter = 10;
    arap_data.with_dynamics = false;

    //precomputation
    if (!igl::arap_precomputation(V, F, 3, b, arap_data)) {
        std::cerr << "precomputation failed\n";
        return 1;
    }
    
    Eigen::MatrixXd res = V;
    //Solve
    if (!igl::arap_solve(bc, arap_data, res)) {
        std::cerr << "ARAP solve failed\n";
        return 1;
    }

    //dispaly/visualization that is copied from main.cpp 
    int n = static_cast<int>(V.rows());
    int m = static_cast<int>(F.rows());

    Eigen::MatrixXd V_original_show = V;
    Eigen::MatrixXd V_deformed_show = res;

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

    points.row(0) = bc.row(0);
    points.row(1) = bc.row(1);

    // shift points to the deformed mesh side
    points.col(0).array() += shift;

    Eigen::MatrixXd point_colors(2, 3);

    // blue = fixed vertex
    point_colors.row(0) = Eigen::RowVector3d(0.0, 0.0, 1.0);

    // red = moved handle vertex
    point_colors.row(1) = Eigen::RowVector3d(1.0, 0.0, 0.0);

    viewer.data().add_points(points, point_colors);

    viewer.launch();
}