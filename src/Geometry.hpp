#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>
#include <map> // 新增

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
    virtual Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const = 0;

    // 新增：计算面上形状函数的积分 int_f phi_i dS
    // 返回一个 map: 节点全局索引 -> 积分值 (weight)
    virtual std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const = 0;
};

// ... (TriFace 和 PolygonFace 的声明同样需要更新) ...
class TriFace : public VEMFace {
public:
    TriFace(const std::vector<int>& nodes) : VEMFace(nodes) {}
    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};

class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes) : VEMFace(nodes) {}
    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};