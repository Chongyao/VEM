#include "VEMMesh.hpp"
#include <fstream>
#include <iostream>
// --- 新增：拷贝构造函数实现 ---
VEMMesh::VEMMesh(const VEMMesh& other) 
    : nodes_(other.nodes_), 
      elements_(other.elements_) 
{
    // 深拷贝 faces
    faces_.reserve(other.faces_.size());
    for (const auto& face_ptr : other.faces_) {
        // 利用 clone() 创建副本
        faces_.push_back(face_ptr->clone());
    }
}

// --- 新增：赋值运算符实现 ---
VEMMesh& VEMMesh::operator=(const VEMMesh& other) {
    if (this != &other) {
        nodes_ = other.nodes_;
        elements_ = other.elements_;
        
        // 清空并重建 faces
        faces_.clear();
        faces_.reserve(other.faces_.size());
        for (const auto& face_ptr : other.faces_) {
            faces_.push_back(face_ptr->clone());
        }
    }
    return *this;
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