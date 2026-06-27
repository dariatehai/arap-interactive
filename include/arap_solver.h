#pragma once

#include "mesh.h"
#include "sparse_builder.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <vector>
#include <unordered_map>

namespace arap {

class ARAPSolver {
public:
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using VertexMatrix = Eigen::MatrixXd;
    using FaceMatrix = Eigen::MatrixXi;

    ARAPSolver() = default;

    // Initialize ARAP from existing Mesh class
    void initialize(const Mesh& mesh);

    // Set fixed / handle vertices.
    // constraint_indices: vertex ids
    // constraint_positions: target positions, size k x 3
    void set_constraints(
        const std::vector<int>& constraint_indices,
        const Eigen::MatrixXd& constraint_positions
    );

    // If only handle positions change but handle vertex ids are the same,
    // this avoids refactorizing the global matrix.
    void update_constraint_positions(
        const Eigen::MatrixXd& constraint_positions
    );

    // Run local-global ARAP iterations
    void solve(int iterations = 5);

    const Eigen::MatrixXd& deformed_vertices() const {
        return V_;
    }

    const Eigen::MatrixXd& original_vertices() const {
        return V0_;
    }

private:
    void precompute_global_solve();

    void local_step();
    void global_step();

    bool is_constraint_vertex(int vertex_id) const;

private:
    // Original and deformed vertex positions
    Eigen::MatrixXd V0_;
    Eigen::MatrixXd V_;

    // Faces are not directly used by ARAP after neighbor precomputation, but useful for exporting a deformed mesh
    Eigen::MatrixXi F_;

    int n_vertices_ = 0;

    // Reuse Mesh::WeightedNeighbor from your Mesh class for local solver
    std::vector<std::vector<Mesh::WeightedNeighbor>> weighted_neighbors_;

    // One local rotation per vertex
    std::vector<Eigen::Matrix3d> rotations_;

    // ARAP Laplacian-like global matrix
    SparseMatrix L_;

    // Constraint data
    std::vector<int> constraint_indices_;
    Eigen::MatrixXd constraint_positions_;

    std::vector<bool> is_constraint_;

    // Free vertex mapping
    std::vector<int> free_indices_; //unknown vertices we need to solve
    std::unordered_map<int, int> global_to_free_;
    std::unordered_map<int, int> global_to_constraint_;

    // Reduced system:
    //
    // L_ff * V_free + L_fc * V_constraint = b_free
    //
    // Therefore:
    //
    // L_ff * V_free = b_free - L_fc * V_constraint
    SparseMatrix L_ff_;
    SparseMatrix L_fc_;

    Eigen::SimplicialLLT<SparseMatrix> solver_;

    bool initialized_ = false;
    bool constraints_ready_ = false;
};

}