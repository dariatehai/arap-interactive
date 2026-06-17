#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <vector>
#include <string>

namespace arap {

/**
 * @brief Mesh representation with vertices, faces, and spatial data
 */
class Mesh {
public:
    using VertexMatrix = Eigen::MatrixXd;
    using FaceMatrix = Eigen::MatrixXi;
    using SparseMatrix = Eigen::SparseMatrix<double>;

    /**
     * @brief Default constructor - creates empty mesh
     */
    Mesh() = default;

    /**
     * @brief Constructor with vertex and face data
     * @param vertices Nx3 matrix of vertex positions
     * @param faces Mx3 matrix of face indices
     */
    Mesh(const VertexMatrix& vertices, const FaceMatrix& faces);

    // Getters
    const VertexMatrix& vertices() const { return vertices_; }
    const FaceMatrix& faces() const { return faces_; }
    int num_vertices() const { return vertices_.rows(); }
    int num_faces() const { return faces_.rows(); }

    // Setters
    void set_vertices(const VertexMatrix& vertices) { vertices_ = vertices; }
    void set_faces(const FaceMatrix& faces) { faces_ = faces; }

    /**
     * @brief Compute per-vertex normals from face data
     * @return Nx3 matrix of vertex normals
     */
    VertexMatrix compute_vertex_normals() const;

    /**
     * @brief Build the cotangent Laplacian matrix
     * @return NxN sparse matrix (cotangent Laplacian)
     */
    SparseMatrix build_cotangent_laplacian() const;

    /**
     * @brief Build the mass matrix (diagonal matrix with vertex areas)
     * @return NxN sparse diagonal matrix
     */
    SparseMatrix build_mass_matrix() const;

    /**
     * @brief Compute vertex areas (using barycentric cell areas)
     * @return N vector of vertex areas
     */
    Eigen::VectorXd compute_vertex_areas() const;

    /**
     * @brief Clear mesh data
     */
    void clear();

    /**
     * @brief Check if mesh is valid (non-empty)
     */
    bool is_valid() const { return vertices_.rows() > 0 && faces_.rows() > 0; }

private:
    VertexMatrix vertices_;  // Nx3 matrix
    FaceMatrix faces_;       // Mx3 matrix
};

}  // namespace arap
