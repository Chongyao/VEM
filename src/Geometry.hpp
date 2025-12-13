#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <map>

struct GeometricProps {
    double area;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
};

class VEMFace {
public:
    std::vector<int> node_indices;

    VEMFace(const std::vector<int>& nodes) : node_indices(nodes) {}
    virtual ~VEMFace() = default;

    // --- 新增：虚拟克隆接口，用于深拷贝 ---
    virtual std::unique_ptr<VEMFace> clone() const = 0;

    virtual GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const = 0;
    virtual Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const = 0;
    virtual std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const = 0;
};

class TriFace : public VEMFace {
public:
    TriFace(const std::vector<int>& nodes) : VEMFace(nodes) {}
    
    // --- 新增实现 ---
    std::unique_ptr<VEMFace> clone() const override;

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};

class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes) : VEMFace(nodes) {}

    // --- 新增实现 ---
    std::unique_ptr<VEMFace> clone() const override;

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};