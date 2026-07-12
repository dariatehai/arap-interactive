#include "ply_parser.h"
#include "arap_solver.h"

#include <igl/opengl/glfw/Viewer.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void make_uv_sphere(
    const Eigen::RowVector3d& center,
    double radius,
    Eigen::MatrixXd& sphere_vertices,
    Eigen::MatrixXi& sphere_faces)
{
    constexpr int latitude_segments = 12;
    constexpr int longitude_segments = 20;
    constexpr double pi = 3.14159265358979323846;

    sphere_vertices.resize((latitude_segments + 1) * longitude_segments, 3);
    for (int latitude = 0; latitude <= latitude_segments; ++latitude) {
        const double phi = pi * latitude / latitude_segments;
        for (int longitude = 0; longitude < longitude_segments; ++longitude) {
            const double theta = 2.0 * pi * longitude / longitude_segments;
            const int index = latitude * longitude_segments + longitude;
            sphere_vertices.row(index) = center + radius * Eigen::RowVector3d(
                std::sin(phi) * std::cos(theta),
                std::cos(phi),
                std::sin(phi) * std::sin(theta));
        }
    }

    sphere_faces.resize(2 * latitude_segments * longitude_segments, 3);
    int face = 0;
    for (int latitude = 0; latitude < latitude_segments; ++latitude) {
        for (int longitude = 0; longitude < longitude_segments; ++longitude) {
            const int next_longitude = (longitude + 1) % longitude_segments;
            const int a = latitude * longitude_segments + longitude;
            const int b = latitude * longitude_segments + next_longitude;
            const int c = (latitude + 1) * longitude_segments + longitude;
            const int d = (latitude + 1) * longitude_segments + next_longitude;
            sphere_faces.row(face++) << a, c, b;
            sphere_faces.row(face++) << b, c, d;
        }
    }
}

} // namespace

int main()
{
    try {
        arap::Mesh mesh = arap::PLYParser::load("data/Simplified_Armadillo.ply");
        if (!mesh.is_valid()) {
            std::cerr << "Error: Loaded mesh is empty." << std::endl;
            return 1;
        }

        Eigen::MatrixXd V0 = mesh.vertices();
        const Eigen::MatrixXi& F = mesh.faces();

        // For now there is exactly one editable vertex. Change this index to
        // choose a different handle; no vertex picking is performed.
        const int handle_vertex = static_cast<int>(V0.rows()) / 16;
        const int fixed_vertex = 0;

        // A single fixed point still permits the mesh to rotate rigidly around
        // it. Anchor its one-ring neighborhood to remove that ambiguity.
        std::vector<int> fixed_vertices = {fixed_vertex};
        for (int face = 0; face < F.rows(); ++face) {
            bool contains_fixed_vertex = false;
            for (int corner = 0; corner < 3; ++corner) {
                contains_fixed_vertex |= F(face, corner) == fixed_vertex;
            }
            if (!contains_fixed_vertex) {
                continue;
            }
            for (int corner = 0; corner < 3; ++corner) {
                const int vertex = F(face, corner);
                if (vertex != handle_vertex &&
                    std::find(fixed_vertices.begin(), fixed_vertices.end(), vertex) ==
                        fixed_vertices.end()) {
                    fixed_vertices.push_back(vertex);
                }
            }
        }

        // Initialize and factorize the existing ARAP solver once. During
        // dragging, only the constraint positions need to be updated.
        arap::ARAPSolver solver;
        solver.initialize(mesh);

        std::vector<int> constraint_indices = fixed_vertices;
        constraint_indices.push_back(handle_vertex);
        const int handle_constraint_row =
            static_cast<int>(constraint_indices.size()) - 1;

        Eigen::MatrixXd constraint_positions(constraint_indices.size(), 3);
        for (int row = 0; row < static_cast<int>(fixed_vertices.size()); ++row) {
            constraint_positions.row(row) = V0.row(fixed_vertices[row]);
        }
        constraint_positions.row(handle_constraint_row) = V0.row(handle_vertex);
        solver.set_constraints(constraint_indices, constraint_positions);

        const double handle_radius =
            0.01 * (V0.colwise().maxCoeff() -
                     V0.colwise().minCoeff()).norm();
        Eigen::MatrixXd handle_vertices;
        Eigen::MatrixXi handle_faces;
        make_uv_sphere(
            V0.row(handle_vertex),
            handle_radius,
            handle_vertices,
            handle_faces);

        // Render the armadillo and handle as one viewer mesh. This prevents
        // libigl from treating the small sphere as a separate camera target.
        Eigen::MatrixXd display_vertices(
            V0.rows() + handle_vertices.rows(), 3);
        display_vertices.topRows(V0.rows()) = V0;
        display_vertices.bottomRows(handle_vertices.rows()) = handle_vertices;

        Eigen::MatrixXi display_faces(F.rows() + handle_faces.rows(), 3);
        display_faces.topRows(F.rows()) = F;
        display_faces.bottomRows(handle_faces.rows()) =
            handle_faces.array() + V0.rows();

        Eigen::MatrixXd display_colors(display_vertices.rows(), 3);
        display_colors.topRows(V0.rows()).rowwise() =
            Eigen::RowVector3d(0.75, 0.75, 0.75);
        display_colors.bottomRows(handle_vertices.rows()).rowwise() =
            Eigen::RowVector3d(1.0, 0.0, 0.0);

        igl::opengl::glfw::Viewer viewer;
        const int mesh_id = viewer.data().id;
        viewer.data().set_mesh(display_vertices, display_faces);
        viewer.data().set_colors(display_colors);
        viewer.data().show_lines = true;

        bool dragging = false;
        float drag_depth = 0.0f;
        Eigen::Vector2f drag_offset = Eigen::Vector2f::Zero();

        // Fit against the original mesh so even the handle's small radius is
        // excluded from the initial camera calculation.
        viewer.callback_init =
            [&](igl::opengl::glfw::Viewer& v) {
                v.core().align_camera_center(V0, F);
                return false;
            };

        viewer.callback_mouse_down =
            [&](igl::opengl::glfw::Viewer& v, int button, int /*modifier*/) {
                if (button != 0) {
                    return false;
                }

                const Eigen::Vector3f handle_position =
                    V0.row(handle_vertex).transpose().cast<float>();
                const Eigen::Vector3f screen_position = igl::project(
                    handle_position,
                    v.core().view,
                    v.core().proj,
                    v.core().viewport);

                const float mouse_screen_y =
                    v.core().viewport(1) + v.core().viewport(3) -
                    v.current_mouse_y;
                const Eigen::Vector2f mouse_position(
                    static_cast<float>(v.current_mouse_x), mouse_screen_y);

                // Consume the click only when it lands on the visible handle.
                // Otherwise libigl receives it and rotates the camera normally.
                constexpr float handle_click_radius_pixels = 24.0f;
                if ((screen_position.head<2>() - mouse_position).norm() >
                    handle_click_radius_pixels) {
                    return false;
                }

                drag_depth = screen_position.z();
                drag_offset = screen_position.head<2>() - mouse_position;
                dragging = true;
                return true;
            };

        viewer.callback_mouse_move =
            [&](igl::opengl::glfw::Viewer& v, int mouse_x, int mouse_y) {
                if (!dragging) {
                    return false;
                }

                // libigl uses a bottom-left origin for screen coordinates,
                // while the mouse callback uses a top-left origin.
                const float screen_y =
                    v.core().viewport(1) + v.core().viewport(3) - mouse_y;
                const Eigen::Vector3f screen_position(
                    static_cast<float>(mouse_x) + drag_offset.x(),
                    screen_y + drag_offset.y(),
                    drag_depth);

                const Eigen::Vector3f world_position = igl::unproject(
                    screen_position,
                    v.core().view,
                    v.core().proj,
                    v.core().viewport);

                // The fixed vertex remains at its original position; the red
                // handle follows the mouse. ARAP computes all other vertices.
                constraint_positions.row(handle_constraint_row) =
                    world_position.cast<double>().transpose();
                solver.update_constraint_positions(constraint_positions);
                solver.solve(50);
                V0 = solver.deformed_vertices();

                make_uv_sphere(
                    constraint_positions.row(handle_constraint_row),
                    handle_radius,
                    handle_vertices,
                    handle_faces);

                display_vertices.topRows(V0.rows()) = V0;
                display_vertices.bottomRows(handle_vertices.rows()) =
                    handle_vertices;
                v.data(mesh_id).set_vertices(display_vertices);
                v.data(mesh_id).compute_normals();
                return true;
            };

        viewer.callback_mouse_up =
            [&](igl::opengl::glfw::Viewer&, int button, int /*modifier*/) {
                if (button == 0 && dragging) {
                    dragging = false;
                    return true;
                }
                return false;
            };

        std::cout << "Hold the left mouse button over the red sphere and drag it.\n"
                  << "Left-drag elsewhere rotates the camera." << std::endl;

        viewer.launch();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
