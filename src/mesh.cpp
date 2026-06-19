#include "mesh.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace arap {

Mesh::Mesh(const VertexMatrix& vertices, const FaceMatrix& faces)
    : vertices_(vertices), faces_(faces) {
    if (vertices.cols() != 3) {
        throw std::invalid_argument("Vertices must be Nx3");
    }
    if (faces.cols() != 3) {
        throw std::invalid_argument("Faces must be Mx3");
    }
}

Mesh::VertexMatrix Mesh::compute_vertex_normals() const {
    VertexMatrix normals = VertexMatrix::Zero(vertices_.rows(), 3);

    // For each face, compute face normal and add to adjacent vertices
    for (int i = 0; i < faces_.rows(); ++i) {
        int v0 = faces_(i, 0);
        int v1 = faces_(i, 1);
        int v2 = faces_(i, 2);

        Eigen::Vector3d e1 = vertices_.row(v1) - vertices_.row(v0);
        Eigen::Vector3d e2 = vertices_.row(v2) - vertices_.row(v0);
        Eigen::Vector3d face_normal = e1.cross(e2);

        normals.row(v0) += face_normal;
        normals.row(v1) += face_normal;
        normals.row(v2) += face_normal;
    }

    // Normalize
    for (int i = 0; i < normals.rows(); ++i) {
        double norm = normals.row(i).norm();
        if (norm > 1e-10) {
            normals.row(i) /= norm;
        }
    }

    return normals;
}

Eigen::VectorXd Mesh::compute_vertex_areas() const {
    Eigen::VectorXd areas = Eigen::VectorXd::Zero(vertices_.rows());

    // Compute barycentric cell areas (1/3 of each triangle area to its vertices)
    for (int i = 0; i < faces_.rows(); ++i) {
        int v0 = faces_(i, 0);
        int v1 = faces_(i, 1);
        int v2 = faces_(i, 2);

        Eigen::Vector3d e1 = vertices_.row(v1) - vertices_.row(v0);
        Eigen::Vector3d e2 = vertices_.row(v2) - vertices_.row(v0);
        double triangle_area = 0.5 * e1.cross(e2).norm();

        // Distribute triangle area equally to its three vertices
        areas(v0) += triangle_area / 3.0;
        areas(v1) += triangle_area / 3.0;
        areas(v2) += triangle_area / 3.0;
    }

    return areas;
}

Mesh::SparseMatrix Mesh::build_cotangent_laplacian() const {
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;

    int n = vertices_.rows();

    // For each edge, compute cotangent weights
    for (int i = 0; i < faces_.rows(); ++i) {
        int v0 = faces_(i, 0);
        int v1 = faces_(i, 1);
        int v2 = faces_(i, 2);

        // Helper lambda to compute cotangent of angle at vertex opposite to edge
        auto compute_cotangent = [this](int v_opposite, int v_edge1, int v_edge2) {
            Eigen::Vector3d u = vertices_.row(v_edge1) - vertices_.row(v_opposite);
            Eigen::Vector3d v = vertices_.row(v_edge2) - vertices_.row(v_opposite);

            double cos_angle = u.dot(v);
            double sin_angle = u.cross(v).norm();

            if (sin_angle < 1e-10) return 0.0;
            return cos_angle / sin_angle;
        };

        // Edge v0-v1, opposite angle at v2
        double cot_v2 = compute_cotangent(v2, v0, v1);
        triplets.push_back(Triplet(v0, v1, cot_v2 / 2.0));
        triplets.push_back(Triplet(v1, v0, cot_v2 / 2.0));

        // Edge v1-v2, opposite angle at v0
        double cot_v0 = compute_cotangent(v0, v1, v2);
        triplets.push_back(Triplet(v1, v2, cot_v0 / 2.0));
        triplets.push_back(Triplet(v2, v1, cot_v0 / 2.0));

        // Edge v2-v0, opposite angle at v1
        double cot_v1 = compute_cotangent(v1, v2, v0);
        triplets.push_back(Triplet(v2, v0, cot_v1 / 2.0));
        triplets.push_back(Triplet(v0, v2, cot_v1 / 2.0));
    }

    // Build matrix and extract diagonal for negative diagonal entries
    SparseMatrix L(n, n);
    L.setFromTriplets(triplets.begin(), triplets.end());

    // Create Laplacian: L = -A (where A is the cotangent matrix)
    // Diagonal entries are negative sum of off-diagonal entries
    for (int i = 0; i < n; ++i) {
        L.coeffRef(i, i) = -L.row(i).sum();
    }

    return L;
}

Mesh::SparseMatrix Mesh::build_mass_matrix() const {
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;

    Eigen::VectorXd areas = compute_vertex_areas();

    // Create diagonal matrix with vertex areas
    for (int i = 0; i < areas.size(); ++i) {
        if (areas(i) > 0) {
            triplets.push_back(Triplet(i, i, areas(i)));
        }
    }

    SparseMatrix M(vertices_.rows(), vertices_.rows());
    M.setFromTriplets(triplets.begin(), triplets.end());

    return M;
}

void Mesh::clear() {
    vertices_.resize(0, 3);
    faces_.resize(0, 3);
    weighted_neighbors_.clear();
}

void Mesh::build_weighted_neighbors()   {
    weighted_neighbors_.clear();
    const int n = vertices_.rows();
    weighted_neighbors_.resize(n);

    if (n == 0 || faces_.rows() == 0) {
        return;
    }
    //Temp edge-weight container: E(i,j) -> w_ij 
    std::map<std::pair<int, int>, double> edge_weights;
    //Create new edge considering the undirected edges E(i,j) = E(j,i)
    auto make_edge_key = [](int a, int b) {
    if (a < b) {
        return std::make_pair(a, b);
    }
        return std::make_pair(b, a);
    };
    // Helper lambda to compute cotangent of angle at vertex opposite to edge
    auto compute_cotangent = [this](int v_opposite, int v_edge1, int v_edge2) {
        Eigen::Vector3d u = vertices_.row(v_edge1) - vertices_.row(v_opposite);
        Eigen::Vector3d v = vertices_.row(v_edge2) - vertices_.row(v_opposite);

        double cos_angle = u.dot(v);
        double sin_angle = u.cross(v).norm();

        if (sin_angle < 1e-10) return 0.0;
        return cos_angle / sin_angle;
    };
    /*
    accumulate the cotangent weight to each edge:
    for boundary only consider one angle;
    if two triganles share an edge, their contributions accumulate
    */
    auto add_edge_weight = [&](int a, int b, double weight) {
        auto key = make_edge_key(a, b);
        edge_weights[key] += weight;
    };
    //iterate over each faces 
    for (int f = 0; f < faces_.rows(); ++f) {
        int v0 = faces_(f, 0);
        int v1 = faces_(f, 1);
        int v2 = faces_(f, 2);

        if (v0 < 0 || v0 >= n || v1 < 0 || v1 >= n || v2 < 0 || v2 >= n) {
            throw std::runtime_error("Face contains invalid vertex index");
        }

        if (v0 == v1 || v1 == v2 || v2 == v0) {
            continue;
        }
        //the weight of an edge uses the angle opposite to that edge
        double cot_v2 = compute_cotangent(v2, v0, v1);//E(v0, v1) is opposite to v2
        double cot_v0 = compute_cotangent(v0, v1, v2);//E(v1, v2) is opposite to v0
        double cot_v1 = compute_cotangent(v1, v2, v0);//E(v2, v0) is opposite to v1

        add_edge_weight(v0, v1, 0.5 * cot_v2);
        add_edge_weight(v1, v2, 0.5 * cot_v0);
        add_edge_weight(v2, v0, 0.5 * cot_v1);
    }
    //convert temp edge-weight map into neighbor list
    for (const auto& entry : edge_weights) {
        int i = entry.first.first;
        int j = entry.first.second;
        double wij = entry.second;

        weighted_neighbors_[i].push_back({j, wij});
        weighted_neighbors_[j].push_back({i, wij});
    }
}   
}  // namespace arap
