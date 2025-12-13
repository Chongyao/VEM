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

// 实现 PolyhedronElement 的几何计算
void PolyhedronElement::updateGeometricProps(const Eigen::MatrixXd& nodes, 
                                             const std::vector<std::unique_ptr<VEMFace>>& faces) {
    volume = 0.0;
    Eigen::Vector3d moment_1st = Eigen::Vector3d::Zero();

    // 1. 计算参考点 X_ref (简单的顶点平均)
    // 这是一个鲁棒性技巧：对于凸多面体，任何点都可以；
    // 对于星形多面体，内核点最好；均值点通常落在内核中。
    Eigen::Vector3d x_ref = Eigen::Vector3d::Zero();
    int node_count = 0;
    // 这里为了简单，我们遍历面的所有点（会有重复，但对均值影响不大，或者可以用 std::set 去重）
    // 为了极致性能，后面可以优化，现在先求面的重心的平均
    for (int fid : face_indices) {
        // 使用面的质心来计算体参考点也是一种稳健策略
        auto props = faces[fid]->computeProps(nodes);
        x_ref += props.centroid;
    }
    if (!face_indices.empty()) x_ref /= face_indices.size();

    // 2. 累加金字塔
    for (int fid : face_indices) {
        // 获取面的属性
        auto face_props = faces[fid]->computeProps(nodes);
        
        // 金字塔高度 (投影到法向)
        // 注意：假设面法向是朝外的 (Outward normal)
        double h = (face_props.centroid - x_ref).dot(face_props.normal);
        
        // 金字塔体积 = 1/3 * BaseArea * Height
        double pyr_vol = (1.0 / 3.0) * face_props.area * h;

        // 金字塔质心 = 3/4 * FaceCentroid + 1/4 * Apex
        Eigen::Vector3d pyr_centroid = 0.75 * face_props.centroid + 0.25 * x_ref;

        // 累加
        volume += pyr_vol;
        moment_1st += pyr_centroid * pyr_vol;
    }

    if (std::abs(volume) > 1e-12) {
        centroid = moment_1st / volume;
    } else {
        centroid = x_ref;
    }
}