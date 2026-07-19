#include "ply_parser.h"
#include "arap_solver.h"

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/project.h>
#include <igl/unproject.h>
#include <Eigen/Core>
#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
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

struct Handle {
    int vertex = -1;
    Eigen::RowVector3d position = Eigen::RowVector3d::Zero();
};

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

        const int fixed_vertex = 0;

        // Pin a small patch instead of one vertex; otherwise the model can
        // still spin around the fixed point.
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
                if (std::find(fixed_vertices.begin(), fixed_vertices.end(), vertex) ==
                        fixed_vertices.end()) {
                    fixed_vertices.push_back(vertex);
                }
            }
        }

        arap::ARAPSolver solver;
        solver.initialize(mesh);

        std::vector<Handle> handles;
        std::vector<int> constraint_indices;
        Eigen::MatrixXd constraint_positions;

        auto update_constraints = [&]() {
            constraint_indices = fixed_vertices;
            constraint_indices.reserve(fixed_vertices.size() + handles.size());
            for (const Handle& handle : handles) {
                constraint_indices.push_back(handle.vertex);
            }

            constraint_positions.resize(
                static_cast<int>(constraint_indices.size()),
                3);
            for (int row = 0; row < static_cast<int>(fixed_vertices.size()); ++row) {
                constraint_positions.row(row) = mesh.vertices().row(fixed_vertices[row]);
            }
            for (int row = 0; row < static_cast<int>(handles.size()); ++row) {
                constraint_positions.row(
                    static_cast<int>(fixed_vertices.size()) + row) =
                    handles[row].position;
            }

            solver.set_constraints(constraint_indices, constraint_positions);
        };

        update_constraints();

        const double handle_radius =
            0.01 * (V0.colwise().maxCoeff() -
                     V0.colwise().minCoeff()).norm();
        Eigen::MatrixXd handle_vertices;
        Eigen::MatrixXi handle_faces;

        auto update_handle_mesh = [&]() {
            constexpr int sphere_vertex_count = (12 + 1) * 20;
            constexpr int sphere_face_count = 2 * 12 * 20;

            handle_vertices.resize(
                static_cast<int>(handles.size()) * sphere_vertex_count,
                3);
            handle_faces.resize(
                static_cast<int>(handles.size()) * sphere_face_count,
                3);

            for (int i = 0; i < static_cast<int>(handles.size()); ++i) {
                Eigen::MatrixXd sphere_vertices;
                Eigen::MatrixXi sphere_faces;
                make_uv_sphere(
                    handles[i].position,
                    handle_radius,
                    sphere_vertices,
                    sphere_faces);

                handle_vertices.middleRows(
                    i * sphere_vertex_count,
                    sphere_vertex_count) = sphere_vertices;
                handle_faces.middleRows(
                    i * sphere_face_count,
                    sphere_face_count) =
                    sphere_faces.array() + i * sphere_vertex_count;
            }
        };

        Eigen::MatrixXd display_vertices;
        Eigen::MatrixXi display_faces;
        Eigen::MatrixXd display_colors;

        auto update_display_mesh = [&]() {
            update_handle_mesh();

            display_vertices.resize(V0.rows() + handle_vertices.rows(), 3);
            display_vertices.topRows(V0.rows()) = V0;
            if (handle_vertices.rows() > 0) {
                display_vertices.bottomRows(handle_vertices.rows()) =
                    handle_vertices;
            }

            display_faces.resize(F.rows() + handle_faces.rows(), 3);
            display_faces.topRows(F.rows()) = F;
            if (handle_faces.rows() > 0) {
                display_faces.bottomRows(handle_faces.rows()) =
                    handle_faces.array() + V0.rows();
            }

            display_colors.resize(display_vertices.rows(), 3);
            display_colors.topRows(V0.rows()).rowwise() =
                Eigen::RowVector3d(0.75, 0.75, 0.75);
            if (handle_vertices.rows() > 0) {
                display_colors.bottomRows(handle_vertices.rows()).rowwise() =
                    Eigen::RowVector3d(1.0, 0.0, 0.0);
            }
        };

        update_display_mesh();

        igl::opengl::glfw::Viewer viewer;
        const int mesh_id = viewer.data().id;
        viewer.data().set_mesh(display_vertices, display_faces);
        viewer.data().set_colors(display_colors);
        viewer.data().show_lines = true;

        auto refresh_viewer_mesh =
            [&](igl::opengl::glfw::Viewer& v, bool update_faces) {
                update_display_mesh();
                if (update_faces) {
                    v.data(mesh_id).clear();
                    v.data(mesh_id).set_mesh(display_vertices, display_faces);
                    v.data(mesh_id).set_colors(display_colors);
                    v.data(mesh_id).show_lines = true;
                } else {
                    v.data(mesh_id).set_vertices(display_vertices);
                    v.data(mesh_id).set_colors(display_colors);
                }
                v.data(mesh_id).compute_normals();
            };

        bool dragging = false;
        int dragged_handle = -1;
        float drag_depth = 0.0f;
        Eigen::Vector2f drag_offset = Eigen::Vector2f::Zero();
        std::string vertex_id_text = "481";
        std::string handle_status = "No handles yet.";

        igl::opengl::glfw::imgui::ImGuiPlugin plugin;
        viewer.plugins.push_back(&plugin);
        igl::opengl::glfw::imgui::ImGuiMenu menu;
        plugin.widgets.push_back(&menu);

        auto push_handle = [&](int vertex_id) {
            if (vertex_id < 0 || vertex_id >= V0.rows()) {
                handle_status = "Vertex id out of range.";
                return;
            }
            if (std::find(fixed_vertices.begin(), fixed_vertices.end(), vertex_id) !=
                fixed_vertices.end()) {
                handle_status = "That vertex is fixed.";
                return;
            }
            const auto duplicate = std::find_if(
                handles.begin(),
                handles.end(),
                [&](const Handle& handle) {
                    return handle.vertex == vertex_id;
            });
            if (duplicate != handles.end()) {
                handle_status = "Already a handle.";
                return;
            }

            handles.push_back({vertex_id, V0.row(vertex_id)});
            update_constraints();
            handle_status = "Added handle " + std::to_string(vertex_id) + ".";
            refresh_viewer_mesh(viewer, true);
        };

        auto pop_handle = [&]() {
            if (handles.empty()) {
                handle_status = "No handle to pop.";
                return;
            }
            const int vertex_id = handles.back().vertex;
            handles.pop_back();
            dragged_handle = -1;
            dragging = false;
            update_constraints();
            solver.solve(5);
            V0 = solver.deformed_vertices();
            handle_status = "Removed handle " + std::to_string(vertex_id) + ".";
            refresh_viewer_mesh(viewer, true);
        };

        auto clear_handles = [&]() {
            handles.clear();
            dragged_handle = -1;
            dragging = false;
            solver.initialize(mesh);
            V0 = mesh.vertices();
            update_constraints();
            handle_status = "Cleared all handles.";
            refresh_viewer_mesh(viewer, true);
        };

        menu.callback_draw_viewer_menu = [&]() {
            menu.draw_viewer_menu();

            if (ImGui::CollapsingHeader("Handles", ImGuiTreeNodeFlags_DefaultOpen)) {
                ImGui::InputText("Vertex ID", vertex_id_text);

                if (ImGui::Button("Push", ImVec2(-1, 0))) {
                    try {
                        size_t parsed = 0;
                        const int vertex_id = std::stoi(vertex_id_text, &parsed);
                        if (parsed != vertex_id_text.size()) {
                            handle_status = "Vertex id must be an integer.";
                        } else {
                            push_handle(vertex_id);
                        }
                    } catch (const std::exception&) {
                        handle_status = "Vertex id must be an integer.";
                    }
                }

                if (ImGui::Button("Pop", ImVec2(-1, 0))) {
                    pop_handle();
                }

                if (ImGui::Button("Clear", ImVec2(-1, 0))) {
                    clear_handles();
                }

                ImGui::Text("Handles: %d", static_cast<int>(handles.size()));
                if (!handles.empty()) {
                    std::string handle_ids = "Stack:";
                    for (const Handle& handle : handles) {
                        handle_ids += " " + std::to_string(handle.vertex);
                    }
                    ImGui::TextWrapped("%s", handle_ids.c_str());
                }
                ImGui::TextWrapped("%s", handle_status.c_str());
            }
        };

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

                const float mouse_screen_y =
                    v.core().viewport(1) + v.core().viewport(3) -
                    v.current_mouse_y;
                const Eigen::Vector2f mouse_position(
                    static_cast<float>(v.current_mouse_x), mouse_screen_y);

                // Leave normal camera controls alone unless the click is on a handle.
                constexpr float handle_click_radius_pixels = 24.0f;
                float best_distance = handle_click_radius_pixels;
                dragged_handle = -1;
                Eigen::Vector3f best_screen_position =
                    Eigen::Vector3f::Zero();
                for (int i = 0; i < static_cast<int>(handles.size()); ++i) {
                    const Eigen::Vector3f handle_position =
                        handles[i].position.transpose().cast<float>();
                    const Eigen::Vector3f screen_position = igl::project(
                        handle_position,
                        v.core().view,
                        v.core().proj,
                        v.core().viewport);
                    const float distance =
                        (screen_position.head<2>() - mouse_position).norm();
                    if (distance <= best_distance) {
                        best_distance = distance;
                        dragged_handle = i;
                        best_screen_position = screen_position;
                    }
                }

                if (dragged_handle < 0) {
                    return false;
                }

                drag_depth = best_screen_position.z();
                drag_offset = best_screen_position.head<2>() - mouse_position;
                dragging = true;
                return true;
            };

        viewer.callback_mouse_move =
            [&](igl::opengl::glfw::Viewer& v, int mouse_x, int mouse_y) {
                if (!dragging) {
                    return false;
                }

                // Mouse y is top-left based; libigl's screen coords are bottom-left based.
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

                handles[dragged_handle].position =
                    world_position.cast<double>().transpose();
                constraint_positions.row(
                    static_cast<int>(fixed_vertices.size()) + dragged_handle) =
                    handles[dragged_handle].position;
                solver.update_constraint_positions(constraint_positions);
                solver.solve(5);
                V0 = solver.deformed_vertices();

                refresh_viewer_mesh(v, false);
                return true;
            };

        viewer.callback_mouse_up =
            [&](igl::opengl::glfw::Viewer&, int button, int /*modifier*/) {
                if (button == 0 && dragging) {
                    dragging = false;
                    dragged_handle = -1;
                    return true;
                }
                return false;
            };

        viewer.launch();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
