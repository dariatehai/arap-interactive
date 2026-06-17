#pragma once

#include "mesh.h"
#include <string>

namespace arap {

/**
 * @brief Parser for PLY (Polygon File Format) files
 */
class PLYParser {
public:
    /**
     * @brief Load a mesh from a PLY file
     * @param filename Path to the PLY file
     * @return Mesh object with loaded data
     * @throws std::runtime_error if file cannot be read or is invalid
     */
    static Mesh load(const std::string& filename);

    /**
     * @brief Save a mesh to a PLY file
     * @param filename Output file path
     * @param mesh Mesh object to save
     * @throws std::runtime_error if file cannot be written
     */
    static void save(const std::string& filename, const Mesh& mesh);

private:
    // Header parsing
    enum class DataType { ASCII, BINARY_LITTLE_ENDIAN, BINARY_BIG_ENDIAN };

    struct Header {
        DataType format;
        int vertex_count = 0;
        int face_count = 0;
        bool has_normals = false;
        bool has_colors = false;
    };

    static Header parse_header(std::ifstream& file);
    static Mesh parse_ascii(std::ifstream& file, const Header& header);
    static Mesh parse_binary(std::ifstream& file, const Header& header, bool little_endian);
};

}  // namespace arap
