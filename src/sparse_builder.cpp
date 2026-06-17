#include "sparse_builder.h"
#include <cmath>

namespace arap {

SparseBuilder::SparseBuilder(int rows, int cols) : rows_(rows), cols_(cols) {}

void SparseBuilder::add(int row, int col, double value) {
    triplets_.push_back(Triplet(row, col, value));
}

void SparseBuilder::set(int row, int col, double value) {
    // Remove duplicates would require more complex logic
    // For now, add() can handle duplicates during build()
    triplets_.push_back(Triplet(row, col, value));
}

SparseBuilder::SparseMatrix SparseBuilder::build() {
    SparseMatrix result(rows_, cols_);
    result.setFromTriplets(triplets_.begin(), triplets_.end());
    return result;
}

void SparseBuilder::clear() {
    triplets_.clear();
}

namespace sparse {

Eigen::SparseMatrix<double> to_sparse(const Eigen::MatrixXd& dense, double threshold) {
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;

    for (int i = 0; i < dense.rows(); ++i) {
        for (int j = 0; j < dense.cols(); ++j) {
            if (std::abs(dense(i, j)) > threshold) {
                triplets.push_back(Triplet(i, j, dense(i, j)));
            }
        }
    }

    Eigen::SparseMatrix<double> result(dense.rows(), dense.cols());
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

Eigen::MatrixXd to_dense(const Eigen::SparseMatrix<double>& sparse) {
    return Eigen::MatrixXd(sparse);
}

double frobenius_norm(const Eigen::SparseMatrix<double>& mat) {
    double norm = 0.0;
    for (int k = 0; k < mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
            norm += it.value() * it.value();
        }
    }
    return std::sqrt(norm);
}

}  // namespace sparse

}  // namespace arap
