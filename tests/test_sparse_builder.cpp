#include <gtest/gtest.h>
#include "sparse_builder.h"
#include <cmath>

using namespace arap;

class SparseBuilderTest : public ::testing::Test {
protected:
    SparseBuilder builder{5, 5};
};

TEST_F(SparseBuilderTest, AddElements) {
    builder.add(0, 0, 1.0);
    builder.add(1, 1, 2.0);
    builder.add(0, 1, 3.0);

    auto mat = builder.build();
    EXPECT_EQ(mat.rows(), 5);
    EXPECT_EQ(mat.cols(), 5);
    EXPECT_EQ(mat.coeff(0, 0), 1.0);
    EXPECT_EQ(mat.coeff(1, 1), 2.0);
    EXPECT_EQ(mat.coeff(0, 1), 3.0);
}

TEST_F(SparseBuilderTest, BuildSparseMatrix) {
    for (int i = 0; i < 5; ++i) {
        builder.add(i, i, double(i + 1));
    }

    auto mat = builder.build();
    EXPECT_EQ(mat.nonZeros(), 5);

    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(mat.coeff(i, i), double(i + 1));
    }
}

TEST_F(SparseBuilderTest, Clear) {
    builder.add(0, 0, 1.0);
    EXPECT_GT(builder.nnz(), 0);

    builder.clear();
    EXPECT_EQ(builder.nnz(), 0);

    auto mat = builder.build();
    EXPECT_EQ(mat.nonZeros(), 0);
}

class SparseUtilTest : public ::testing::Test {};

TEST_F(SparseUtilTest, DenseToSparse) {
    Eigen::MatrixXd dense(3, 3);
    dense << 1.0, 0.0, 2.0,
             0.0, 3.0, 0.0,
             4.0, 0.0, 5.0;

    auto sparse = sparse::to_sparse(dense);
    EXPECT_EQ(sparse.rows(), 3);
    EXPECT_EQ(sparse.cols(), 3);
    EXPECT_EQ(sparse.coeff(0, 0), 1.0);
    EXPECT_EQ(sparse.coeff(0, 2), 2.0);
    EXPECT_EQ(sparse.coeff(1, 1), 3.0);
}

TEST_F(SparseUtilTest, SparseToEnse) {
    Eigen::MatrixXd original(3, 3);
    original << 1.0, 0.0, 2.0,
                0.0, 3.0, 0.0,
                4.0, 0.0, 5.0;

    auto sparse = sparse::to_sparse(original);
    auto dense = sparse::to_dense(sparse);

    EXPECT_TRUE(dense.isApprox(original));
}

TEST_F(SparseUtilTest, FrobeniusNorm) {
    SparseBuilder builder(3, 3);
    builder.add(0, 0, 3.0);
    builder.add(1, 1, 4.0);

    auto mat = builder.build();
    double norm = sparse::frobenius_norm(mat);

    // sqrt(3^2 + 4^2) = sqrt(9 + 16) = sqrt(25) = 5
    EXPECT_NEAR(norm, 5.0, 1e-10);
}

TEST_F(SparseUtilTest, FrobeniusNormZero) {
    SparseBuilder builder(3, 3);
    auto mat = builder.build();
    double norm = sparse::frobenius_norm(mat);

    EXPECT_EQ(norm, 0.0);
}
