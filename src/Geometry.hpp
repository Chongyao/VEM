#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>

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

    virtual GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const = 0;

    // 修改：返回 Vector4d [Area, Mx, My, Mz]
    // 这是计算体积分和 VEM 投影算子的基础
    virtual Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const = 0;
};

class TriFace : public VEMFace {
public:
    TriFace(const std::vector<int>& nodes) : VEMFace(nodes) {}
    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
};

class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes) : VEMFace(nodes) {}
    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
};