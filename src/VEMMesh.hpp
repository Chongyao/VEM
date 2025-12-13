#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include "Geometry.hpp"

// 单元定义：仅存储面的索引
struct PolyhedronElement {
    int id;
    std::vector<int> face_indices; // 构成本单元的面的 ID 列表
    
    // 缓存数据（可选，VEM 计算时填充）
    double volume;
    Eigen::Vector3d centroid;
};

class VEMMesh {
private:
    // SoA 数据：3行 N列
    Eigen::MatrixXd nodes_; 
    
    // 拓扑存储
    std::vector<std::unique_ptr<VEMFace>> faces_;
    std::vector<PolyhedronElement> elements_;

public:
    // 初始化节点 (3 x N 矩阵)
    void setNodes(const Eigen::MatrixXd& V) {
        nodes_ = V;
    }

    // 添加三角形面
    int addTriFace(const std::vector<int>& ids) {
        faces_.push_back(std::make_unique<TriFace>(ids));
        return faces_.size() - 1;
    }

    // 添加多边形面 (通用)
    int addPolygonFace(const std::vector<int>& ids) {
        faces_.push_back(std::make_unique<PolygonFace>(ids));
        return faces_.size() - 1;
    }

    // 添加单元 (通过面的索引)
    void addElement(const std::vector<int>& face_ids) {
        PolyhedronElement elem;
        elem.id = elements_.size();
        elem.face_indices = face_ids;
        // 在这里可以立即计算体积，或者lazy compute
        elements_.push_back(elem);
    }

    const Eigen::MatrixXd& getNodes() const { return nodes_; }
    const VEMFace* getFace(int id) const { return faces_[id].get(); }

    // ==========================================
    // 轻量级 VTK (Unstructured Grid .vtu) 输出
    // ==========================================
    void exportToVTK(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out) { std::cerr << "Error opening file " << filename << std::endl; return; }

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
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
        // VTK_POLYHEDRON (Type 42) 格式极其特殊：
        // Format: [num_faces, num_pts_face0, id... id, num_pts_face1, id... id, ...]
        
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<int> types;
        
        int current_offset = 0;

        for (const auto& elem : elements_) {
            types.push_back(42); // VTK_POLYHEDRON
            
            // 统计该 poly 的数据大小
            // header = 1 (num_faces)
            // for each face: 1 (num_pts) + num_pts
            int cell_data_size = 1; 
            for (int fid : elem.face_indices) {
                cell_data_size += 1 + faces_[fid]->node_indices.size();
            }
            
            offsets.push_back(current_offset + cell_data_size);
            current_offset += cell_data_size;

            // 填充 Connectivity 数据
            connectivity.push_back(elem.face_indices.size()); // Number of faces
            for (int fid : elem.face_indices) {
                const auto& node_ids = faces_[fid]->node_indices;
                connectivity.push_back(node_ids.size());      // Number of points in face
                for (int nid : node_ids) {
                    connectivity.push_back(nid);              // Point IDs
                }
            }
        }

        out << "      <Cells>\n";
        
        // Connectivity
        out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (int v : connectivity) out << v << " ";
        out << "\n        </DataArray>\n";

        // Offsets
        out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (int v : offsets) out << v << " ";
        out << "\n        </DataArray>\n";

        // Types
        out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (int v : types) out << v << " ";
        out << "\n        </DataArray>\n";
        
        // Face Stream (VTK Polyhedron 需要额外的 faces 数据，但在 XML 格式中，
        // 实际上 connectivity 承载了这一信息。上述 connectivity 写法是针对 XML Polyhedron 的标准写法)
        // 注意：这是 VTK 5.0+ 的 XML 格式处理 Polyhedron 的方式。

        out << "      </Cells>\n";
        out << "    </Piece>\n";
        out << "  </UnstructuredGrid>\n";
        out << "</VTKFile>\n";
        
        std::cout << "Exported VTK to " << filename << std::endl;
    }

};