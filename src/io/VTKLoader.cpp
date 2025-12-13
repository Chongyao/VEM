#include "VTKLoader.hpp"
#include <tinyxml2.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

using namespace tinyxml2;

namespace vem::io {

// 辅助函数：将字符串 "1.0 2.0 3.0 ..." 解析为 vector
template <typename T>
std::vector<T> parseArray(const char* text) {
    std::vector<T> res;
    if (!text) return res;
    std::stringstream ss(text);
    T val;
    while (ss >> val) res.push_back(val);
    return res;
}

std::unique_ptr<VEMMesh> VTKLoader::load(const std::string& filename) const {
    XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
        std::cerr << "[VTKLoader] Error: Failed to load file " << filename << std::endl;
        return nullptr;
    }

    XMLElement* root = doc.FirstChildElement("VTKFile");
    if (!root) return nullptr;
    
    XMLElement* grid = root->FirstChildElement("UnstructuredGrid");
    if (!grid) return nullptr;
    
    XMLElement* piece = grid->FirstChildElement("Piece");
    if (!piece) return nullptr;

    auto mesh = std::make_unique<VEMMesh>();

    // 1. 解析节点 (Points)
    XMLElement* points = piece->FirstChildElement("Points");
    if (points) {
        XMLElement* da = points->FirstChildElement("DataArray");
        if (da && da->GetText()) {
            std::vector<double> raw = parseArray<double>(da->GetText());
            // 假设是 3D
            int n_points = raw.size() / 3;
            Eigen::MatrixXd V(3, n_points);
            for (int i = 0; i < n_points; ++i) {
                V(0, i) = raw[3*i+0];
                V(1, i) = raw[3*i+1];
                V(2, i) = raw[3*i+2];
            }
            mesh->setNodes(V);
        }
    }

    // 2. 解析单元 (Cells)
    XMLElement* cells = piece->FirstChildElement("Cells");
    if (cells) {
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<int> types;
        std::vector<int> faces;
        std::vector<int> faceoffsets;

        XMLElement* da = cells->FirstChildElement("DataArray");
        while (da) {
            const char* name = da->Attribute("Name");
            if (name) {
                std::string s(name);
                if (s == "connectivity") connectivity = parseArray<int>(da->GetText());
                else if (s == "offsets") offsets = parseArray<int>(da->GetText());
                else if (s == "types") types = parseArray<int>(da->GetText());
                else if (s == "faces") faces = parseArray<int>(da->GetText());
                else if (s == "faceoffsets") faceoffsets = parseArray<int>(da->GetText());
            }
            da = da->NextSiblingElement("DataArray");
        }

        // 构建 VEM 单元
        int n_elems = types.size();
        int conn_idx = 0; // connectivity 指针

        for (int i = 0; i < n_elems; ++i) {
            int type = types[i];
            
            // 计算当前单元在 connectivity 中的节点数
            // 如果是第一个单元，offset 就是 size
            // 否则 size = offset[i] - offset[i-1]
            int end_offset = offsets[i];
            int start_offset = (i == 0) ? 0 : offsets[i-1];
            int n_cell_nodes = end_offset - start_offset;

            // 获取当前单元的节点 ID 列表 (仅用于标准单元)
            std::vector<int> node_ids;
            for(int k=0; k<n_cell_nodes; ++k) {
                node_ids.push_back(connectivity[start_offset + k]);
            }

            // 处理不同类型
            if (type == 10) { // VTK_TETRA
                // 4 nodes: 0,1,2,3
                // Faces (VTK standard): (0,1,2), (0,2,3), (0,3,1), (1,3,2) -> wait, check outward normal
                // VTK Tet order: n0, n1, n2, n3.
                // Face 0: 0-1-2 (n3 on top? normal = (1-0)x(2-0). if n3 above, normal points down(in)? check)
                // 通常 VEM 不依赖严格的面顺序，只要闭合即可，几何引擎会自动计算体积符号。
                // 简单起见，添加4个面
                std::vector<int> f_ids;
                f_ids.push_back(mesh->addTriFace({node_ids[0], node_ids[2], node_ids[1]})); // Flip for safety?
                f_ids.push_back(mesh->addTriFace({node_ids[0], node_ids[1], node_ids[3]}));
                f_ids.push_back(mesh->addTriFace({node_ids[0], node_ids[3], node_ids[2]}));
                f_ids.push_back(mesh->addTriFace({node_ids[1], node_ids[2], node_ids[3]}));
                mesh->addElement(f_ids);
            }
            else if (type == 12) { // VTK_HEXAHEDRON
                // 8 nodes
                // Faces: 6 quads
                // 0-3-2-1, 4-5-6-7, 0-1-5-4, 1-2-6-5, 2-3-7-6, 3-0-4-7 (VTK ordering)
                std::vector<int> f_ids;
                f_ids.push_back(mesh->addPolygonFace({node_ids[0], node_ids[3], node_ids[2], node_ids[1]})); // Bottom
                f_ids.push_back(mesh->addPolygonFace({node_ids[4], node_ids[5], node_ids[6], node_ids[7]})); // Top
                f_ids.push_back(mesh->addPolygonFace({node_ids[0], node_ids[1], node_ids[5], node_ids[4]})); // Front
                f_ids.push_back(mesh->addPolygonFace({node_ids[1], node_ids[2], node_ids[6], node_ids[5]})); // Right
                f_ids.push_back(mesh->addPolygonFace({node_ids[2], node_ids[3], node_ids[7], node_ids[6]})); // Back
                f_ids.push_back(mesh->addPolygonFace({node_ids[3], node_ids[0], node_ids[4], node_ids[7]})); // Left
                mesh->addElement(f_ids);
            }
            else if (type == 42) { // VTK_POLYHEDRON
                // 需要查看 faces 数组
                // faceoffsets[i] 指向 faces 数组的起始位置
                if (faceoffsets.empty()) {
                    std::cerr << "Error: VTK_POLYHEDRON found but no 'faces' or 'faceoffsets' array." << std::endl;
                    continue;
                }
                
                int f_start = (i==0) ? 0 : faceoffsets[i-1];
                // int f_end = faceoffsets[i]; // not needed strictly if data is consistent
                
                // faces array structure: [nCellFaces, nFace0Nodes, id..., nFace1Nodes, id...]
                int n_cell_faces = faces[f_start];
                int ptr = f_start + 1;
                
                std::vector<int> elem_face_ids;
                for(int f=0; f<n_cell_faces; ++f) {
                    int n_face_nodes = faces[ptr++];
                    std::vector<int> face_nodes;
                    for(int k=0; k<n_face_nodes; ++k) {
                        face_nodes.push_back(faces[ptr++]);
                    }
                    elem_face_ids.push_back(mesh->addPolygonFace(face_nodes));
                }
                mesh->addElement(elem_face_ids);
            }
        }
    }

    std::cout << "[VTKLoader] Loaded " << mesh->getNumNodes() << " nodes, " 
              << mesh->getNumElements() << " elements." << std::endl;
    return mesh;
}

} // namespace vem::io