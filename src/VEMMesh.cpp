#include "VEMMesh.hpp"
#include <fstream>
#include <iostream>
// ==========================================
// 轻量级 VTK (Unstructured Grid .vtu) 输出 (修复版)
// ==========================================
void VEMMesh::exportToVTK(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) { std::cerr << "Error opening file " << filename << std::endl; return; }

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        out << "  <UnstructuredGrid>\n";
        out << "    <Piece NumberOfPoints=\"" << nodes_.cols() << "\" NumberOfCells=\"" << elements_.size() << "\">\n";

        // 1. 输出节点坐标
        out << "      <Points>\n";
        out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int i = 0; i < nodes_.cols(); ++i) {
            out << nodes_(0, i) << " " << nodes_(1, i) << " " << nodes_(2, i) << " ";
        }
        out << "\n        </DataArray>\n";
        out << "      </Points>\n";

        // 2. 输出单元 (Cells)
        // XML 格式下，VTK_POLYHEDRON 必须包含 faces 和 faceoffsets 数组
        
        std::vector<int> connectivity;  // 仅存储构成单元的点 ID
        std::vector<int> offsets;       // connectivity 的偏移量
        std::vector<int> types;         // 单元类型 (42)
        std::vector<int> faces;         // [nFaces, nPtsF1, p1.., nPtsF2, p2..]
        std::vector<int> faceoffsets;   // faces 的偏移量
        
        int current_offset = 0;
        int current_face_offset = 0;

        for (const auto& elem : elements_) {
            types.push_back(42); // VTK_POLYHEDRON
            
            // --- 构建 faces 和 faceoffsets (你原来的 connectivity 逻辑应该放在这里) ---
            int cell_face_data_size = 1; // 头部包含 1 个数：面的数量
            faces.push_back(elem.face_indices.size()); // Number of faces

            // 临时收集该单元所有涉及的点，用于填充 connectivity
            // (虽然 faces 里已经有了点，但 connectivity 依然需要列出属于该单元的点，供 Paraview 索引)
            std::vector<int> cell_node_ids;

            for (int fid : elem.face_indices) {
                const auto& node_ids = faces_[fid]->node_indices;
                
                // 填充 faces 数组
                faces.push_back(node_ids.size());      // Number of points in face
                for (int nid : node_ids) {
                    faces.push_back(nid);              // Point IDs
                    cell_node_ids.push_back(nid);      // 收集点ID用于 connectivity
                }
                
                cell_face_data_size += 1 + node_ids.size();
            }
            
            current_face_offset += cell_face_data_size;
            faceoffsets.push_back(current_face_offset);

            // --- 构建 connectivity 和 offsets ---
            // 注意：connectivity 只需要列出点 ID，不需要“面数”或“点数”这些头信息
            for(int nid : cell_node_ids) {
                connectivity.push_back(nid);
            }
            current_offset += cell_node_ids.size();
            offsets.push_back(current_offset);
        }

        out << "      <Cells>\n";
        
        // 数组 1: Connectivity (只含点 ID)
        out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (int v : connectivity) out << v << " ";
        out << "\n        </DataArray>\n";

        // 数组 2: Offsets (对应 Connectivity)
        out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (int v : offsets) out << v << " ";
        out << "\n        </DataArray>\n";

        // 数组 3: Types
        out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (int v : types) out << v << " ";
        out << "\n        </DataArray>\n";

        // 数组 4: Faces (关键！你原来的代码把这个写进了 connectivity)
        out << "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">\n";
        for (int v : faces) out << v << " ";
        out << "\n        </DataArray>\n";

        // 数组 5: FaceOffsets (关键！)
        out << "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">\n";
        for (int v : faceoffsets) out << v << " ";
        out << "\n        </DataArray>\n";

        out << "      </Cells>\n";
        out << "    </Piece>\n";
        out << "  </UnstructuredGrid>\n";
        out << "</VTKFile>\n";
        
        std::cout << "Exported VTK to " << filename << std::endl;
    }