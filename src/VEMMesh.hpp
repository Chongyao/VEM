#pragma once
#include <vector>
#include <string>

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
    void exportToVTK(const std::string& filename) const;

};