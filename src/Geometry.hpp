#pragma once
#include <vector>
#include <Eigen/Dense>
#include <memory>

// 定义几何属性结构体
struct GeometricProps {
    double area;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
    // VEM 需要计算的投影矩阵 D 相关的积分可以在这里扩展
};

// 抽象面基类
class VEMFace {
public:
    std::vector<int> node_indices; // 逆时针顺序

    VEMFace(const std::vector<int>& nodes) : node_indices(nodes) {}
    virtual ~VEMFace() = default;

    // 核心几何计算：传入节点坐标矩阵 (3 x N)，利用 SoA 高效计算
    virtual GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const = 0;

    // 计算 VEM 投影所需的单项式积分 (int_F phi * n ds)
    // 这里先简化，只留接口
    virtual void integrateMonomials(const Eigen::MatrixXd& all_nodes, 
                                    Eigen::MatrixXd& integrals) const = 0;
};

// 三角形面实现
class TriFace : public VEMFace {
public:
    TriFace(const std::vector<int>& nodes) : VEMFace(nodes) {}

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;

    void integrateMonomials(const Eigen::MatrixXd& all_nodes, Eigen::MatrixXd& integrals) const override;
};

// 多边形面（用于 Cut-VEM 和四边形）
class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes) : VEMFace(nodes) {}

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    
    void integrateMonomials(const Eigen::MatrixXd& all_nodes, Eigen::MatrixXd& integrals) const override {
        // 后续实现
    }
};