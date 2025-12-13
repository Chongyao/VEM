#pragma once
#include <vector>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include "Geometry.hpp"

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

public:
    void setNodes(const Eigen::MatrixXd& V) { nodes_ = V; }

    int addTriFace(const std::vector<int>& ids) {
        faces_.push_back(std::make_unique<TriFace>(ids));
        return faces_.size() - 1;
    }

    int addPolygonFace(const std::vector<int>& ids) {
        faces_.push_back(std::make_unique<PolygonFace>(ids));
        return faces_.size() - 1;
    }

    void addElement(const std::vector<int>& face_ids) {
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

    void exportToVTK(const std::string& filename) const;
};