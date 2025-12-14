#include "VTKWriter.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace vem::io {

// 辅助：写入网格几何与拓扑的公共部分
// 由于要写 XML 结构，直接由 save 函数控制比较好。
// 这里我们实现带数据的版本，如果不带数据，就把 u 传空或者分开写。

void VTKWriter::save(const VEMMesh& mesh, const std::string& filename) const
{
    // 调用带数据的版本，传入空向量表示无数据
    save(mesh, Eigen::VectorXd(), filename);
}

void VTKWriter::save(const VEMMesh& mesh, const Eigen::VectorXd& u, const std::string& filename) const
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "[VTKWriter] Error: Could not open " << filename << "\n";
        return;
    }

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << mesh.getNumNodes() << "\" NumberOfCells=\"" << mesh.getNumElements() << "\">\n";

    // --- Points ---
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    const auto& nodes = mesh.getNodes();
    for (int i = 0; i < mesh.getNumNodes(); ++i) {
        out << nodes(0, i) << " " << nodes(1, i) << " " << nodes(2, i) << "\n";
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // --- Cells ---
    out << "      <Cells>\n";

    // 1. connectivity
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        const auto& el = mesh.getElement(i);

        // VTK Polyhedron connectivity format is implicit in standard VTK XML,
        // but for general polyhedra we usually use Face stream.
        // HOWEVER, to keep it simple and compatible with typical readers:
        // Standard cells (Tet, Hex) act normally.
        // Polyhedron (Type 42) uses 'faces' array. The 'connectivity' array for Polyhedron
        // usually lists the point IDs but strictly speaking for Type 42,
        // the geometry is defined in the 'faces' array.
        // ParaView expects 'connectivity' to still contain points?
        // Actually for VTK_POLYHEDRON, 'connectivity' is ignored by some parsers if 'faces' is present,
        // but valid VTK requires it to list all unique points of the cell.
        // Let's just list all nodes of all faces.

        // Collect unique nodes for this element
        std::vector<int> unique_nodes;
        for (int fid : el.face_indices) {
            const auto* face = mesh.getFace(fid);
            for (int nid : face->node_indices)
                unique_nodes.push_back(nid);
        }
        // Deduplicate
        std::sort(unique_nodes.begin(), unique_nodes.end());
        unique_nodes.erase(std::unique(unique_nodes.begin(), unique_nodes.end()), unique_nodes.end());

        for (int nid : unique_nodes)
            out << nid << " ";
        out << "\n";
    }
    out << "        </DataArray>\n";

    // 2. offsets
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int current_offset = 0;
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        const auto& el = mesh.getElement(i);

        // Re-calculate unique nodes count
        std::vector<int> unique_nodes;
        for (int fid : el.face_indices) {
            const auto* face = mesh.getFace(fid);
            for (int nid : face->node_indices)
                unique_nodes.push_back(nid);
        }
        std::sort(unique_nodes.begin(), unique_nodes.end());
        unique_nodes.erase(std::unique(unique_nodes.begin(), unique_nodes.end()), unique_nodes.end());

        current_offset += unique_nodes.size();
        out << current_offset << "\n";
    }
    out << "        </DataArray>\n";

    // 3. types
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        out << "42\n"; // VTK_POLYHEDRON
    }
    out << "        </DataArray>\n";

    // 4. faces (Required for VTK_POLYHEDRON)
    // Format: [nFaces, nFace0Pts, id..., nFace1Pts, id..., ...]
    out << "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">\n";
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        const auto& el = mesh.getElement(i);

        out << el.face_indices.size() << " "; // Number of faces
        for (int fid : el.face_indices) {
            const auto* face = mesh.getFace(fid);
            out << face->node_indices.size() << " ";
            for (int nid : face->node_indices)
                out << nid << " ";
        }
        out << "\n";
    }
    out << "        </DataArray>\n";

    // 5. faceoffsets (Required for VTK_POLYHEDRON)
    out << "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">\n";
    int current_face_offset = 0;
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        const auto& el = mesh.getElement(i);

        // Count total ints for this cell in the 'faces' array
        int cell_ints = 1; // The face count itself
        for (int fid : el.face_indices) {
            const auto* face = mesh.getFace(fid);
            cell_ints += 1 + face->node_indices.size(); // 1 for numPts + pts
        }
        current_face_offset += cell_ints;
        out << current_face_offset << "\n";
    }
    out << "        </DataArray>\n";

    out << "      </Cells>\n";

    // --- Point Data (Displacement) ---
    if (u.size() == 3 * mesh.getNumNodes()) {
        out << "      <PointData Vectors=\"Displacement\">\n";
        out << "        <DataArray type=\"Float64\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < mesh.getNumNodes(); ++i) {
            out << u(3 * i) << " " << u(3 * i + 1) << " " << u(3 * i + 2) << "\n";
        }
        out << "        </DataArray>\n";
        out << "      </PointData>\n";
    }

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    out.close();
}

} // namespace vem::io
