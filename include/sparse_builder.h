#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <utility>

namespace arap {

/**
 * @brief Utility class for building sparse matrices efficiently
 */
class SparseBuilder {
public:
    using Triplet = Eigen::Triplet<double>;
    using SparseMatrix = Eigen::SparseMatrix<double>;

    /**
     * @brief Constructor
     * @param rows Number of rows in the matrix
     * @param cols Number of columns in the matrix
     */
    explicit SparseBuilder(int rows, int cols);

    /**
     * @brief Add a value to the matrix at position (row, col)
     * @param row Row index
     * @param col Column index
     * @param value Value to add
     */
    void add(int row, int col, double value);

    /**
     * @brief Set a value in the matrix at position (row, col) (overwrites)
     * @param row Row index
     * @param col Column index
     * @param value Value to set
     */
    void set(int row, int col, double value);

    /**
     * @brief Build and return the sparse matrix
     * @return Constructed sparse matrix
     */
    SparseMatrix build();

    /**
     * @brief Clear all entries
     */
    void clear();

    /**
     * @brief Get number of non-zero entries
     */
    int nnz() const { return triplets_.size(); }

private:
    int rows_, cols_;
    std::vector<Triplet> triplets_;
};

/**
 * @brief Utility functions for sparse matrix operations
 */
namespace sparse {

/**
 * @brief Convert dense matrix to sparse format
 * @param dense Input dense matrix
     * @param threshold Drop elements below this threshold
     * @return Sparse matrix representation
     */
    Eigen::SparseMatrix<double> to_sparse(
        const Eigen::MatrixXd& dense,
        double threshold = 1e-12);

    /**
     * @brief Convert sparse matrix to dense format
     * @param sparse Input sparse matrix
     * @return Dense matrix representation
     */
    Eigen::MatrixXd to_dense(const Eigen::SparseMatrix<double>& sparse);

    /**
     * @brief Compute Frobenius norm of sparse matrix
     * @param mat Input sparse matrix
     * @return Frobenius norm
     */
    double frobenius_norm(const Eigen::SparseMatrix<double>& mat);

}  // namespace sparse

}  // namespace arap
