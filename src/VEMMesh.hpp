#pragma once
#include "Geometry.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

struct PolyhedronElement {
    int id;
    std::vector<int> face_indices;

    // 缓存数据
    double volume = 0.0;
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();

    // 新增：计算几何属性
    void updateGeometricProps(const Eigen::MatrixXd& nodes,
        const std::vector<std::unique_ptr<VEMFace>>& faces);
};

class VEMMesh {
private:
    Eigen::MatrixXd nodes_;
    std::vector<std::unique_ptr<VEMFace>> faces_;
    std::vector<PolyhedronElement> elements_;

    std::vector<std::vector<int>> face_to_elements_;

public:
    VEMMesh() = default;

    VEMMesh(const VEMMesh& other); // Copy Constructor
    VEMMesh& operator=(const VEMMesh& other); // Copy Assignment

    VEMMesh(VEMMesh&&) = default;
    VEMMesh& operator=(VEMMesh&&) = default;
    void setNodes(const Eigen::MatrixXd& V) { nodes_ = V; }

    int getNumElements() const { return elements_.size(); }

    int addTriFace(const std::vector<int>& ids)
    {
        faces_.push_back(std::make_unique<TriFace>(ids));
        return faces_.size() - 1;
    }

    int addPolygonFace(const std::vector<int>& ids)
    {
        faces_.push_back(std::make_unique<PolygonFace>(ids));
        return faces_.size() - 1;
    }

    void addElement(const std::vector<int>& face_ids)
    {
        PolyhedronElement elem;
        elem.id = elements_.size();
        elem.face_indices = face_ids;
        // 创建时立即计算几何属性
        elem.updateGeometricProps(nodes_, faces_);
        elements_.push_back(elem);
    }

    const Eigen::MatrixXd& getNodes() const { return nodes_; }
    const VEMFace* getFace(int id) const { return faces_[id].get(); }
    const PolyhedronElement& getElement(int id) const { return elements_[id]; }
    int getNumNodes() const { return nodes_.cols(); }

    // [新增] 构建/更新拓扑邻接关系
    // 每次网格发生变动（如加载后、合并后）都需要调用
    void buildTopology();

    // [新增] 获取共享指定面的单元列表
    const std::vector<int>& getElementsOnFace(int face_id) const;

    // [新增] 获取某单元的所有相邻单元 ID
    // 用于 Agglomeration 寻找合并候选者
    std::vector<int> getElementNeighbors(int elem_id) const;

    // [新增] 检查网格是否为流形 (Manifold)
    // 即检查是否所有内部面都恰好被 2 个单元共享
    bool checkManifold() const;
    int getNumFaces() const { return faces_.size(); }
};
