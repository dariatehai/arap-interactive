#include "ply_parser.h"

#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>

#include <iostream>
#include <stdexcept>

int main()
{
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

    return 0;
}