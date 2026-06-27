#include "arap_solver.h"
#include "mesh.h"

#include <gtest/gtest.h>
#include <Eigen/Core>

TEST(ARAPSolverTest, BasicSolveTetrahedron) {
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

    arap::Mesh mesh(vertices, faces);

    arap::ARAPSolver solver;
    solver.initialize(mesh);

    std::vector<int> constraints = {0, 1};

    Eigen::MatrixXd constraint_positions(2, 3);
    constraint_positions.row(0) = vertices.row(0);
    constraint_positions.row(1) = vertices.row(1) + Eigen::RowVector3d(0.2, 0.0, 0.0);

    solver.set_constraints(constraints, constraint_positions);
    solver.solve(5);

    Eigen::MatrixXd result = solver.deformed_vertices();

    EXPECT_EQ(result.rows(), 4);
    EXPECT_EQ(result.cols(), 3);
    EXPECT_TRUE(result.allFinite());

    EXPECT_NEAR((result.row(0) - constraint_positions.row(0)).norm(), 0.0, 1e-8);
    EXPECT_NEAR((result.row(1) - constraint_positions.row(1)).norm(), 0.0, 1e-8);
}