#include "arap_solver.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace arap {

void ARAPSolver::initialize(const Mesh& mesh) {
    V0_ = mesh.vertices();  //original vertices
    V_ = V0_;               //deformed vertices
    F_ = mesh.faces();

    if (V0_.cols() != 3) {
        throw std::runtime_error("ARAPSolver expects vertices to be n x 3.");
    }

    n_vertices_ = static_cast<int>(V0_.rows());

    
    // neighbor lists and corresponding weights
    weighted_neighbors_ = mesh.build_weighted_neighbors();

    if (static_cast<int>(weighted_neighbors_.size()) != n_vertices_) {
        throw std::runtime_error("Weighted neighbor size does not match vertex count.");
    }
    //initialize the R for each vertices
    rotations_.assign(n_vertices_, Eigen::Matrix3d::Identity());
    //create Laplacian-like matrix L for global part
    L_ = -mesh.build_cotangent_laplacian();

    initialized_ = true;
    constraints_ready_ = false;
}

//change constraint vertex IDs and positions
void ARAPSolver::set_constraints(
    const std::vector<int>& constraint_indices,
    const Eigen::MatrixXd& constraint_positions
) {
    if (!initialized_) {
        throw std::runtime_error("Call initialize() before set_constraints().");
    }
    //global solver requires constraints
    if (constraint_indices.empty()) {
        throw std::runtime_error("ARAP needs at least one constraint vertex.");
    }

    if (constraint_positions.rows() != static_cast<int>(constraint_indices.size()) ||
        constraint_positions.cols() != 3) {
        throw std::runtime_error("Constraint positions must be k x 3.");
    }
    //store constraints in the member variable
    constraint_indices_ = constraint_indices;   //indices of handle vertices
    //Small matrix, containing only constraint vertices
    constraint_positions_ = constraint_positions;   

    is_constraint_.assign(n_vertices_, false);
    //unordered_map that maps global vertex index to row index in the [constraint_positions_]
    global_to_constraint_.clear();

    for (int i = 0; i < static_cast<int>(constraint_indices_.size()); ++i) {
        int vertex_id = constraint_indices_[i];

        if (vertex_id < 0 || vertex_id >= n_vertices_) {
            throw std::runtime_error("Constraint vertex index out of range.");
        }

        is_constraint_[vertex_id] = true;
        global_to_constraint_[vertex_id] = i;

        // Initialize constrained vertex position
        V_.row(vertex_id) = constraint_positions_.row(i);
    }

    precompute_global_solve();

    constraints_ready_ = true;
}
//only change the target positions
void ARAPSolver::update_constraint_positions(
    const Eigen::MatrixXd& constraint_positions
) {
    if (!constraints_ready_) {
        throw std::runtime_error("Call set_constraints() before update_constraint_positions().");
    }

    if (constraint_positions.rows() != static_cast<int>(constraint_indices_.size()) ||
        constraint_positions.cols() != 3) {
        throw std::runtime_error("Constraint positions must be k x 3.");
    }

    constraint_positions_ = constraint_positions;

    for (int c = 0; c < static_cast<int>(constraint_indices_.size()); ++c) {
        int vertex_id = constraint_indices_[c];
        V_.row(vertex_id) = constraint_positions_.row(c);
    }
}

void ARAPSolver::precompute_global_solve() {
    free_indices_.clear();
    global_to_free_.clear();

    for (int i = 0; i < n_vertices_; ++i) {
        if (!is_constraint_[i]) {
            int local_id = static_cast<int>(free_indices_.size());
            free_indices_.push_back(i);
            global_to_free_[i] = local_id;
        }
    }

    const int n_free = static_cast<int>(free_indices_.size());
    const int n_constraints = static_cast<int>(constraint_indices_.size());

    SparseBuilder builder_ff(n_free, n_free);
    SparseBuilder builder_fc(n_free, n_constraints);

    /*
     Iterate over global sparse matrix L.
     For each non-zero entry L(row, col), decide whether it belongs to:
     L_ff: free-free block
     L_fc: free-constraint block
    */
    for (int outer = 0; outer < L_.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator it(L_, outer); it; ++it) {
            int row = static_cast<int>(it.row());
            int col = static_cast<int>(it.col());
            double value = it.value();
            
            if (is_constraint_[row]) {
                continue;
            }

            int local_row = global_to_free_.at(row);

            if (is_constraint_[col]) {
                int local_col = global_to_constraint_.at(col);
                builder_fc.add(local_row, local_col, value);
            } else {
                int local_col = global_to_free_.at(col);
                builder_ff.add(local_row, local_col, value);
            }
        }
    }

    L_ff_ = builder_ff.build();
    L_fc_ = builder_fc.build();

    //performs the Cholesky factorization
    solver_.compute(L_ff_);

    if (solver_.info() != Eigen::Success) {
        throw std::runtime_error("Failed to factorize ARAP global matrix.");
    }
}

void ARAPSolver::local_step() {
    //iterate all vertices to update the changes
    for (int i = 0; i < n_vertices_; ++i) {
        //S_i from Equation 5 in the paper
        Eigen::Matrix3d S = Eigen::Matrix3d::Zero();

        for (const auto& neighbor : weighted_neighbors_[i]) {
            int j = neighbor.vertex;
            double w = neighbor.weight;

            Eigen::Vector3d p_ij = V0_.row(i) - V0_.row(j);
            Eigen::Vector3d q_ij = V_.row(i) - V_.row(j);

            // covariance matrix
            S += w * p_ij * q_ij.transpose();
        }

        Eigen::JacobiSVD<Eigen::Matrix3d> svd(
            S,
            Eigen::ComputeFullU | Eigen::ComputeFullV
        );

        Eigen::Matrix3d U = svd.matrixU();
        Eigen::Matrix3d V = svd.matrixV();

        Eigen::Matrix3d R = V * U.transpose();

        // Avoid reflection
        if (R.determinant() < 0.0) {
            V.col(2) *= -1.0;
            R = V * U.transpose();
        }

        rotations_[i] = R;
    }
}

void ARAPSolver::global_step() {
    //Equation 8 in the paper
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n_vertices_, 3);
    //construct the b
    for (int i = 0; i < n_vertices_; ++i) {
        Eigen::Vector3d rhs_i = Eigen::Vector3d::Zero();

        for (const auto& neighbor : weighted_neighbors_[i]) {
            int j = neighbor.vertex;
            double w = neighbor.weight;

            Eigen::Vector3d p_ij = V0_.row(i) - V0_.row(j);

            rhs_i += 0.5 * w * (rotations_[i] + rotations_[j]) * p_ij;
        }

        b.row(i) = rhs_i.transpose();
    }

    const int n_free = static_cast<int>(free_indices_.size());

    //instead of solving all p', we solve only free vertices(unkown vertices)
    Eigen::MatrixXd b_free(n_free, 3);

    for (int local_i = 0; local_i < n_free; ++local_i) {
        int global_i = free_indices_[local_i];
        b_free.row(local_i) = b.row(global_i);
    }

    // L_ff * V_free + L_fc * V_constraint = b_free
    // Therefore:
    // L_ff * V_free = b_free - L_fc * V_constraint

    Eigen::MatrixXd rhs = b_free - L_fc_ * constraint_positions_;
    
    //uses the factorization to solve the linear system
    Eigen::MatrixXd solved_free = solver_.solve(rhs);

    if (solver_.info() != Eigen::Success) {
        throw std::runtime_error("Failed to solve ARAP global step.");
    }

    //update the deformed vertices based on the result
    for (int local_i = 0; local_i < n_free; ++local_i) {
        int global_i = free_indices_[local_i];
        V_.row(global_i) = solved_free.row(local_i);
    }

    // Keep constraints exactly fixed
    for (int c = 0; c < static_cast<int>(constraint_indices_.size()); ++c) {
        int global_i = constraint_indices_[c];
        V_.row(global_i) = constraint_positions_.row(c);
    }
}

void ARAPSolver::solve(int iterations) {
    if (!initialized_) {
        throw std::runtime_error("Call initialize() before solve().");
    }

    if (!constraints_ready_) {
        throw std::runtime_error("Call set_constraints() before solve().");
    }

    for (int iter = 0; iter < iterations; ++iter) {
        local_step();
        global_step();
    }
}

bool ARAPSolver::is_constraint_vertex(int vertex_id) const {
    if (vertex_id < 0 || vertex_id >= n_vertices_) {
        return false;
    }

    return is_constraint_[vertex_id];
}

} // namespace arap