#include <gtest/gtest.h>
#include "mesh.h"
#include <cmath>

using namespace arap;

class MeshTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple tetrahedron
        Eigen::MatrixXd vertices(4, 3);
        vertices << 0, 0, 0,
                    1, 0, 0,
                    0, 1, 0,
                    0, 0, 1;

        Eigen::MatrixXi faces(4, 3);
        faces << 0, 1, 2,
                 0, 1, 3,
                 0, 2, 3,
                 1, 2, 3;

        mesh = Mesh(vertices, faces);
    }

    Mesh mesh;
};

TEST_F(MeshTest, CreationAndGetters) {
    EXPECT_EQ(mesh.num_vertices(), 4);
    EXPECT_EQ(mesh.num_faces(), 4);
    EXPECT_TRUE(mesh.is_valid());
}

TEST_F(MeshTest, VertexNormals) {
    auto normals = mesh.compute_vertex_normals();
    EXPECT_EQ(normals.rows(), 4);
    EXPECT_EQ(normals.cols(), 3);

    // Check that normals are unit vectors
    for (int i = 0; i < normals.rows(); ++i) {
        double norm = normals.row(i).norm();
        EXPECT_GT(norm, 0.0);  // Should be normalized (or zero if no adjacent faces)
    }
}

TEST_F(MeshTest, VertexAreas) {
    auto areas = mesh.compute_vertex_areas();
    EXPECT_EQ(areas.size(), 4);

    // All areas should be positive
    for (int i = 0; i < areas.size(); ++i) {
        EXPECT_GE(areas(i), 0.0);
    }

    // Sum of areas should be positive
    EXPECT_GT(areas.sum(), 0.0);
}

TEST_F(MeshTest, CotangentLaplacian) {
    auto L = mesh.build_cotangent_laplacian();
    EXPECT_EQ(L.rows(), 4);
    EXPECT_EQ(L.cols(), 4);

    // Check that diagonal entries are negative sums of off-diagonal entries
    for (int i = 0; i < 4; ++i) {
        double row_sum = 0;
        for (int j = 0; j < 4; ++j) {
            row_sum += L.coeff(i, j);
        }
        // Row sum should be approximately zero (sum of all entries in row)
        EXPECT_LT(std::abs(row_sum), 1e-10);
    }
}

TEST_F(MeshTest, MassMatrix) {
    auto M = mesh.build_mass_matrix();
    EXPECT_EQ(M.rows(), 4);
    EXPECT_EQ(M.cols(), 4);

    // Mass matrix should be diagonal
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (i != j) {
                EXPECT_EQ(M.coeff(i, j), 0.0);
            }
        }
    }

    // Diagonal entries should be positive
    for (int i = 0; i < 4; ++i) {
        EXPECT_GT(M.coeff(i, i), 0.0);
    }
}

TEST_F(MeshTest, ClearMesh) {
    mesh.clear();
    EXPECT_EQ(mesh.num_vertices(), 0);
    EXPECT_EQ(mesh.num_faces(), 0);
    EXPECT_FALSE(mesh.is_valid());
}

TEST(MeshConstructor, InvalidDimensions) {
    // Vertices with wrong dimensions should throw
    Eigen::MatrixXd bad_vertices(4, 4);  // 4x4 instead of Nx3
    Eigen::MatrixXi faces(0, 3);

    EXPECT_THROW(Mesh(bad_vertices, faces), std::invalid_argument);
}
