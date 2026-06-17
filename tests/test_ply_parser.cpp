#include <gtest/gtest.h>
#include "ply_parser.h"
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace arap;

class PLYParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a simple test PLY file
        test_filename = "test_mesh.ply";
        std::ofstream file(test_filename);
        file << "ply\n";
        file << "format ascii 1.0\n";
        file << "element vertex 4\n";
        file << "property float x\n";
        file << "property float y\n";
        file << "property float z\n";
        file << "element face 2\n";
        file << "property list uchar int vertex_indices\n";
        file << "end_header\n";
        file << "0 0 0\n";
        file << "1 0 0\n";
        file << "0 1 0\n";
        file << "0 0 1\n";
        file << "3 0 1 2\n";
        file << "3 0 1 3\n";
        file.close();
    }

    void TearDown() override {
        std::remove(test_filename.c_str());
    }

    std::string test_filename;
};

TEST_F(PLYParserTest, LoadAsciiPLY) {
    Mesh mesh = PLYParser::load(test_filename);

    EXPECT_EQ(mesh.num_vertices(), 4);
    EXPECT_EQ(mesh.num_faces(), 2);

    // Check vertex positions
    EXPECT_NEAR(mesh.vertices()(0, 0), 0.0, 1e-6);
    EXPECT_NEAR(mesh.vertices()(1, 0), 1.0, 1e-6);
    EXPECT_NEAR(mesh.vertices()(2, 1), 1.0, 1e-6);
    EXPECT_NEAR(mesh.vertices()(3, 2), 1.0, 1e-6);

    // Check face indices
    EXPECT_EQ(mesh.faces()(0, 0), 0);
    EXPECT_EQ(mesh.faces()(0, 1), 1);
    EXPECT_EQ(mesh.faces()(0, 2), 2);
    EXPECT_EQ(mesh.faces()(1, 0), 0);
    EXPECT_EQ(mesh.faces()(1, 1), 1);
    EXPECT_EQ(mesh.faces()(1, 2), 3);
}

TEST_F(PLYParserTest, SaveAsciiPLY) {
    std::string output_filename = "test_output.ply";

    // Create a simple mesh
    Eigen::MatrixXd vertices(3, 3);
    vertices << 0, 0, 0,
                1, 0, 0,
                0, 1, 0;

    Eigen::MatrixXi faces(1, 3);
    faces << 0, 1, 2;

    Mesh mesh(vertices, faces);
    PLYParser::save(output_filename, mesh);

    // Reload and verify
    Mesh loaded_mesh = PLYParser::load(output_filename);
    EXPECT_EQ(loaded_mesh.num_vertices(), 3);
    EXPECT_EQ(loaded_mesh.num_faces(), 1);

    // Clean up
    std::remove(output_filename.c_str());
}

TEST_F(PLYParserTest, InvalidFile) {
    EXPECT_THROW(PLYParser::load("nonexistent_file.ply"), std::runtime_error);
}

TEST_F(PLYParserTest, InvalidPLYFormat) {
    std::string bad_filename = "bad_format.ply";
    std::ofstream file(bad_filename);
    file << "not_ply\n";
    file.close();

    EXPECT_THROW(PLYParser::load(bad_filename), std::runtime_error);

    std::remove(bad_filename.c_str());
}
