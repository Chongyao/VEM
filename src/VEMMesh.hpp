#pragma once
#include "Geometry.hpp" // 必须包含这个，因为它定义了 PolyhedronElement
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <string>
#include <vector>

class VEMMesh {
private:
    Eigen::MatrixXd nodes_;
    std::vector<std::unique_ptr<VEMFace>> faces_;
    // PolyhedronElement 现在定义在 Geometry.hpp 中
    std::vector<PolyhedronElement> elements_;

    // 拓扑邻接信息
    std::vector<std::vector<int>> face_to_elements_;

    // 面去重查找表
    std::map<std::vector<int>, int> face_map_;

public:
    VEMMesh() = default;

    // [Fix] 恢复拷贝构造函数和赋值运算符（深拷贝）
    VEMMesh(const VEMMesh& other);
    VEMMesh& operator=(const VEMMesh& other);

    // 允许移动
    VEMMesh(VEMMesh&&) = default;
    VEMMesh& operator=(VEMMesh&&) = default;

    void setNodes(const Eigen::MatrixXd& V) { nodes_ = V; }

    int getNumElements() const { return elements_.size(); }
    int getNumFaces() const { return faces_.size(); }
    int getNumNodes() const { return nodes_.cols(); }

    const Eigen::MatrixXd& getNodes() const { return nodes_; }
    const VEMFace* getFace(int id) const { return faces_[id].get(); }
    const PolyhedronElement& getElement(int id) const { return elements_[id]; }

    // --- Mesh Construction & Editing ---

    int addTriFace(const std::vector<int>& ids);
    int addPolygonFace(const std::vector<int>& ids);

    int addElement(const std::vector<int>& face_ids);
    void deactivateElement(int elem_id);
    bool isElementActive(int elem_id) const;

    void garbageCollect();
    // --- Topology ---

    void buildTopology();
    const std::vector<int>& getElementsOnFace(int face_id) const;
    std::vector<int> getElementNeighbors(int elem_id) const;
    bool checkManifold() const;
};
