#include "ply_parser.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace arap {

PLYParser::Header PLYParser::parse_header(std::ifstream& file) {
    Header header;
    std::string line;

    // Read magic number
    std::getline(file, line);
    if (line != "ply") {
        throw std::runtime_error("Invalid PLY file: missing 'ply' magic number");
    }

    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "format") {
            std::string format_type;
            iss >> format_type;
            if (format_type == "ascii") {
                header.format = DataType::ASCII;
            } else if (format_type == "binary_little_endian") {
                header.format = DataType::BINARY_LITTLE_ENDIAN;
            } else if (format_type == "binary_big_endian") {
                header.format = DataType::BINARY_BIG_ENDIAN;
            }
        } else if (token == "element") {
            std::string element_type;
            int count;
            iss >> element_type >> count;
            if (element_type == "vertex") {
                header.vertex_count = count;
            } else if (element_type == "face") {
                header.face_count = count;
            }
        } else if (token == "property") {
            // We could parse properties more carefully, but for now we skip them
        } else if (token == "end_header") {
            break;
        }
    }

    if (header.vertex_count == 0) {
        throw std::runtime_error("Invalid PLY file: no vertices specified");
    }

    return header;
}

Mesh PLYParser::parse_ascii(std::ifstream& file, const Header& header) {
    Mesh::VertexMatrix vertices(header.vertex_count, 3);
    Mesh::FaceMatrix faces(header.face_count, 3);

    // Read vertices
    for (int i = 0; i < header.vertex_count; ++i) {
        double x, y, z;
        if (!(file >> x >> y >> z)) {
            throw std::runtime_error("Error reading vertex data");
        }
        vertices(i, 0) = x;
        vertices(i, 1) = y;
        vertices(i, 2) = z;

        // Skip any additional properties (normals, colors, etc.)
        std::string dummy;
        std::getline(file, dummy);
    }

    // Read faces
    for (int i = 0; i < header.face_count; ++i) {
        int num_verts;
        int v0, v1, v2;
        if (!(file >> num_verts >> v0 >> v1 >> v2)) {
            throw std::runtime_error("Error reading face data");
        }
        if (num_verts != 3) {
            throw std::runtime_error("Only triangular faces are supported");
        }
        faces(i, 0) = v0;
        faces(i, 1) = v1;
        faces(i, 2) = v2;
    }

    return Mesh(vertices, faces);
}

Mesh PLYParser::parse_binary(std::ifstream& file, const Header& header, bool little_endian) {
    // Binary parsing would require byte-level reading and endianness handling
    // For now, throw an error - binary support can be added later
    throw std::runtime_error("Binary PLY format not yet supported");
}

Mesh PLYParser::load(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // First pass: read header in text mode
    file.close();
    file.open(filename);

    Header header = parse_header(file);

    if (header.vertex_count == 0 || header.face_count == 0) {
        throw std::runtime_error("Invalid mesh: no vertices or faces");
    }

    switch (header.format) {
        case DataType::ASCII:
            return parse_ascii(file, header);
        case DataType::BINARY_LITTLE_ENDIAN:
            return parse_binary(file, header, true);
        case DataType::BINARY_BIG_ENDIAN:
            return parse_binary(file, header, false);
        default:
            throw std::runtime_error("Unknown PLY format");
    }
}

void PLYParser::save(const std::string& filename, const Mesh& mesh) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }

    const auto& vertices = mesh.vertices();
    const auto& faces = mesh.faces();

    // Write header
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "element vertex " << vertices.rows() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "element face " << faces.rows() << "\n";
    file << "property list uchar int vertex_indices\n";
    file << "end_header\n";

    // Write vertices
    for (int i = 0; i < vertices.rows(); ++i) {
        file << vertices(i, 0) << " " << vertices(i, 1) << " " << vertices(i, 2) << "\n";
    }

    // Write faces
    for (int i = 0; i < faces.rows(); ++i) {
        file << "3 " << faces(i, 0) << " " << faces(i, 1) << " " << faces(i, 2) << "\n";
    }

    file.close();
}

}  // namespace arap
