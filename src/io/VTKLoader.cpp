#include "VTKLoader.hpp"
#include <algorithm> // for transform
#include <iostream>
#include <sstream>
#include <string>
#include <tinyxml2.h>
#include <vector>

using namespace tinyxml2;

namespace vem::io {

// 辅助：转小写
std::string toLower(const std::string& str)
{
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

// 辅助：解析数组
template <typename T>
std::vector<T> parseArray(const char* text)
{
    std::vector<T> res;
    if (!text)
        return res;

    std::string s(text);
    if (s.find_first_not_of(" \t\n\r") == std::string::npos)
        return res;

    std::stringstream ss(s);
    T val;
    while (ss >> val)
        res.push_back(val);
    return res;
}

std::unique_ptr<VEMMesh> VTKLoader::load(const std::string& filename) const
{
    std::cout << "[VTKLoader] Loading " << filename << "..." << std::endl;
    XMLDocument doc;
    if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
        std::cerr << "[VTKLoader] Error: Failed to open file/XML format error." << std::endl;
        return nullptr;
    }

    XMLElement* root = doc.FirstChildElement("VTKFile");
    if (!root)
        return nullptr;

    XMLElement* grid = root->FirstChildElement("UnstructuredGrid");
    if (!grid)
        return nullptr;

    XMLElement* piece = grid->FirstChildElement("Piece");
    if (!piece)
        return nullptr;

    auto mesh = std::make_unique<VEMMesh>();

    // 1. 解析节点 (Points)
    XMLElement* points = piece->FirstChildElement("Points");
    if (points) {
        XMLElement* da = points->FirstChildElement("DataArray");
        if (da) {
            const char* fmt = da->Attribute("format");
            if (fmt && (std::string(fmt) == "binary" || std::string(fmt) == "appended")) {
                std::cerr << "[VTKLoader] CRITICAL ERROR: Binary VTK not supported. Use ASCII." << std::endl;
                return nullptr;
            }
            if (da->GetText()) {
                std::vector<double> raw = parseArray<double>(da->GetText());
                int n_points = raw.size() / 3;
                Eigen::MatrixXd V(3, n_points);
                for (int i = 0; i < n_points; ++i) {
                    V(0, i) = raw[3 * i + 0];
                    V(1, i) = raw[3 * i + 1];
                    V(2, i) = raw[3 * i + 2];
                }
                mesh->setNodes(V);
            }
        }
    }

    // 2. 解析单元 (Cells)
    XMLElement* cells = piece->FirstChildElement("Cells");
    if (cells) {
        // 标准/旧版 Polyhedron 格式
        std::vector<int> connectivity;
        std::vector<int> offsets;
        std::vector<int> types;
        std::vector<int> faces; // Legacy format
        std::vector<int> faceoffsets; // Legacy format

        // VTK 9.0+ 新版 Polyhedron 格式
        std::vector<int> face_connectivity;
        std::vector<int> face_offsets;
        std::vector<int> polyhedron_to_faces;
        std::vector<int> polyhedron_offsets;

        XMLElement* da = cells->FirstChildElement("DataArray");
        while (da) {
            const char* name = da->Attribute("Name");
            if (name) {
                std::string s = toLower(name);

                // Debug: 打印找到的数组名，方便调试
                // std::cout << "Found DataArray: " << s << std::endl;

                if (s == "connectivity")
                    connectivity = parseArray<int>(da->GetText());
                else if (s == "offsets")
                    offsets = parseArray<int>(da->GetText());
                else if (s == "types")
                    types = parseArray<int>(da->GetText());
                // 旧版兼容
                else if (s == "faces")
                    faces = parseArray<int>(da->GetText());
                else if (s == "faceoffsets")
                    faceoffsets = parseArray<int>(da->GetText());
                // 新版兼容 (注意下划线)
                else if (s == "face_connectivity")
                    face_connectivity = parseArray<int>(da->GetText());
                else if (s == "face_offsets")
                    face_offsets = parseArray<int>(da->GetText());
                else if (s == "polyhedron_to_faces")
                    polyhedron_to_faces = parseArray<int>(da->GetText());
                else if (s == "polyhedron_offsets")
                    polyhedron_offsets = parseArray<int>(da->GetText());
            }
            da = da->NextSiblingElement("DataArray");
        }

        // 检查我们处于哪种模式
        bool use_legacy_poly = (!faces.empty() && !faceoffsets.empty());
        bool use_new_poly = (!face_connectivity.empty() && !polyhedron_offsets.empty());

        if (!use_legacy_poly && !use_new_poly) {
            // 如果存在 Polyhedron 类型但没有任何面数据，这是一个问题
            // 但我们先继续，可能只有 Tet 单元
        }

        int n_elems = types.size();
        for (int i = 0; i < n_elems; ++i) {
            int type = types[i];

            // 计算标准单元 (Tet/Hex) 的节点
            int end_offset = (i < offsets.size()) ? offsets[i] : 0;
            int start_offset = (i == 0) ? 0 : offsets[i - 1];
            int n_cell_nodes = end_offset - start_offset;

            std::vector<int> node_ids;
            if (n_cell_nodes > 0 && start_offset < connectivity.size()) {
                for (int k = 0; k < n_cell_nodes; ++k) {
                    node_ids.push_back(connectivity[start_offset + k]);
                }
            }

            if (type == 10) { // VTK_TETRA
                // 即使是 Tet，也可以添加为 VEM 元素 (由4个三角形面组成)
                // 你的代码之前是把 Tet 拆成 4 个 TriFace
                std::vector<int> f_ids;
                f_ids.push_back(mesh->addTriFace({ node_ids[0], node_ids[2], node_ids[1] }));
                f_ids.push_back(mesh->addTriFace({ node_ids[0], node_ids[1], node_ids[3] }));
                f_ids.push_back(mesh->addTriFace({ node_ids[0], node_ids[3], node_ids[2] }));
                f_ids.push_back(mesh->addTriFace({ node_ids[1], node_ids[2], node_ids[3] }));
                mesh->addElement(f_ids);
            } else if (type == 12) { // VTK_HEXAHEDRON
                std::vector<int> f_ids;
                f_ids.push_back(mesh->addPolygonFace({ node_ids[0], node_ids[3], node_ids[2], node_ids[1] }));
                f_ids.push_back(mesh->addPolygonFace({ node_ids[4], node_ids[5], node_ids[6], node_ids[7] }));
                f_ids.push_back(mesh->addPolygonFace({ node_ids[0], node_ids[1], node_ids[5], node_ids[4] }));
                f_ids.push_back(mesh->addPolygonFace({ node_ids[1], node_ids[2], node_ids[6], node_ids[5] }));
                f_ids.push_back(mesh->addPolygonFace({ node_ids[2], node_ids[3], node_ids[7], node_ids[6] }));
                f_ids.push_back(mesh->addPolygonFace({ node_ids[3], node_ids[0], node_ids[4], node_ids[7] }));
                mesh->addElement(f_ids);
            } else if (type == 42) { // VTK_POLYHEDRON
                std::vector<int> elem_face_ids;

                if (use_legacy_poly) {
                    // --- 旧版逻辑 (faces stream) ---
                    // faces array: [nFaces, nPts, id..., nPts, id...]
                    // faceoffsets 指向 faces 数组的起始

                    if (i >= faceoffsets.size())
                        continue;
                    int f_start = (i == 0) ? 0 : faceoffsets[i - 1];
                    if (f_start >= faces.size())
                        continue;

                    int n_cell_faces = faces[f_start];
                    int ptr = f_start + 1;

                    for (int f = 0; f < n_cell_faces; ++f) {
                        int n_face_nodes = faces[ptr++];
                        std::vector<int> face_nodes;
                        for (int k = 0; k < n_face_nodes; ++k)
                            face_nodes.push_back(faces[ptr++]);
                        elem_face_ids.push_back(mesh->addPolygonFace(face_nodes));
                    }
                } else if (use_new_poly) {
                    // --- 新版逻辑 (Decomposed arrays) --- [关键修复]
                    // 1. 确定该单元包含哪些面 (polyhedron_to_faces)
                    int poly_end = (i < polyhedron_offsets.size()) ? polyhedron_offsets[i] : 0;
                    int poly_start = (i == 0) ? 0 : polyhedron_offsets[i - 1];

                    // 遍历该单元的每个面索引
                    for (int p_idx = poly_start; p_idx < poly_end; ++p_idx) {
                        int face_idx = polyhedron_to_faces[p_idx]; // 获取面的全局索引

                        // 2. 查找该面的节点 (face_connectivity)
                        // 注意：face_idx 是相对于 face_offsets 的
                        if (face_idx >= face_offsets.size())
                            continue;

                        int face_end = face_offsets[face_idx];
                        int face_start = (face_idx == 0) ? 0 : face_offsets[face_idx - 1];

                        std::vector<int> face_nodes;
                        for (int k = face_start; k < face_end; ++k) {
                            face_nodes.push_back(face_connectivity[k]);
                        }

                        elem_face_ids.push_back(mesh->addPolygonFace(face_nodes));
                    }
                } else {
                    std::cerr << "[VTKLoader] Error: Cell " << i << " is Polyhedron but no face data found.\n";
                    continue;
                }

                mesh->addElement(elem_face_ids);
            }
        }
    }

    std::cout << "[VTKLoader] Loaded " << mesh->getNumNodes() << " nodes, "
              << mesh->getNumElements() << " elements.\n";
    std::cout << "[VTKLoader] Building mesh topology...\n";
    mesh->buildTopology(); // <--- 新增这行
    mesh->checkManifold(); // <--- 可选：打印一下拓扑状态

    std::cout << "[VTKLoader] Loaded " << mesh->getNumNodes() << " nodes, "
              << mesh->getNumElements() << " elements.\n";
    return mesh;
}

} // namespace vem::io
