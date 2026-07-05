#include <igl/arap.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <iostream>

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

int main() {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //read mesh
    if(!igl::read_triangle_mesh("data/Simplified_Armadillo.ply", V, F)){
        std::cerr << "Failed to read mesh.\n";
        return 1;
    }

    // For now there is exactly one editable vertex. Change this index to
    // choose a different handle; no vertex picking is performed.
    const int handle_vertex = static_cast<int>(V.rows()) / 16;
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

    //number of constraints (size of fixed constraint + 1 editable handle)
    int num_constraints = static_cast<int>(fixed_vertices.size()) + 1;
    //Define constraints b(n): n constraint indices
    Eigen::VectorXi b(num_constraints);
    //Define target position as a n*3 matrix
    Eigen::MatrixXd bc(num_constraints, 3);
    //Add all fixed constraints.
    //same meaning as the line 100-110 in main.cpp
    for (int row = 0; row < static_cast<int>(fixed_vertices.size()); ++row) {
        const int vertex_index = fixed_vertices[row];
        b(row) = vertex_index;
        // Fixed vertices stay at their original positions.
        bc.row(row) = V.row(vertex_index);
    }
    // Put the movable handle in the final constraint row.
    const int handle_constraint_row =  static_cast<int>(fixed_vertices.size());
    b(handle_constraint_row) = handle_vertex;
    bc.row(handle_constraint_row) = V.row(handle_vertex);

    const double handle_radius =
            0.01 * (V.colwise().maxCoeff() -
                     V.colwise().minCoeff()).norm();
        Eigen::MatrixXd handle_vertices;
        Eigen::MatrixXi handle_faces;
        make_uv_sphere(
            V.row(handle_vertex),
            handle_radius,
            handle_vertices,
            handle_faces);
    
    //ARAP refering to Github libigl/include/igl/arap.h
    igl::ARAPData arap_data;
    //configure arap
    arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
    //iteration time is same as in the main.cpp
    arap_data.max_iter = 5;
    arap_data.with_dynamics = false;

    //precomputation
    if (!igl::arap_precomputation(V, F, 3, b, arap_data)) {
        std::cerr << "precomputation failed\n";
        return 1;
    }

    //copied from main.cpp (line 123 - end)
    // Render the armadillo and handle as one viewer mesh. This prevents
    // libigl from treating the small sphere as a separate camera target.
    Eigen::MatrixXd display_vertices(V.rows() + handle_vertices.rows(), 3);
    display_vertices.topRows(V.rows()) = V;
    display_vertices.bottomRows(handle_vertices.rows()) = handle_vertices;

    Eigen::MatrixXi display_faces(F.rows() + handle_faces.rows(), 3);
    display_faces.topRows(F.rows()) = F;
    display_faces.bottomRows(handle_faces.rows()) = handle_faces.array() + V.rows();

    Eigen::MatrixXd display_colors(display_vertices.rows(), 3);
    display_colors.topRows(V.rows()).rowwise() = Eigen::RowVector3d(0.75, 0.75, 0.75);
    display_colors.bottomRows(handle_vertices.rows()).rowwise() = Eigen::RowVector3d(1.0, 0.0, 0.0);
    
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
            v.core().align_camera_center(V, F);
            return false;
        };

    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& v, int button, int /*modifier*/) {
            if (button != 0) {
                return false;
            }

            const Eigen::Vector3f handle_position =
                V.row(handle_vertex).transpose().cast<float>();
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
                bc.row(handle_constraint_row) =
                    world_position.cast<double>().transpose();
                    
                //Solve, the deformed vertices are stored in [res]
                if (!igl::arap_solve(bc, arap_data, V)) {
                    std::cerr << "ARAP solve failed\n";
                    return false;
                }

                make_uv_sphere(
                    bc.row(handle_constraint_row),
                    handle_radius,
                    handle_vertices,
                    handle_faces);

                display_vertices.topRows(V.rows()) = V;
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

        return 0;
}