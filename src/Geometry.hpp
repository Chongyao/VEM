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

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override {
        // 从 SoA 矩阵中获取列向量
        Eigen::Vector3d p0 = all_nodes.col(node_indices[0]);
        Eigen::Vector3d p1 = all_nodes.col(node_indices[1]);
        Eigen::Vector3d p2 = all_nodes.col(node_indices[2]);

        Eigen::Vector3d edge1 = p1 - p0;
        Eigen::Vector3d edge2 = p2 - p0;
        Eigen::Vector3d cross = edge1.cross(edge2);

        GeometricProps props;
        props.area = 0.5 * cross.norm();
        props.normal = cross.normalized();
        props.centroid = (p0 + p1 + p2) / 3.0;
        return props;
    }

    void integrateMonomials(const Eigen::MatrixXd& all_nodes, Eigen::MatrixXd& integrals) const override {
        // 后续实现：三角形上的线性积分
    }
};

// 多边形面（用于 Cut-VEM 和四边形）
class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes) : VEMFace(nodes) {}

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override {
        // 使用 "Shoelace formula" 的 3D 推广
        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        Eigen::Vector3d normal_weighted = Eigen::Vector3d::Zero();
        Eigen::Vector3d p0 = all_nodes.col(node_indices[0]);

        // 简单的扇形剖分求和（假设凸或星形，对于一般多边形需更稳健算法）
        for (size_t i = 1; i < node_indices.size() - 1; ++i) {
            Eigen::Vector3d p1 = all_nodes.col(node_indices[i]);
            Eigen::Vector3d p2 = all_nodes.col(node_indices[i+1]);
            Eigen::Vector3d cross = (p1 - p0).cross(p2 - p0);
            normal_weighted += cross;
            
            // 累加三角形重心 * 面积
            double tri_area = 0.5 * cross.norm();
            Eigen::Vector3d tri_center = (p0 + p1 + p2) / 3.0;
            centroid += tri_center * tri_area; 
        }

        GeometricProps props;
        props.area = 0.5 * normal_weighted.norm();
        props.normal = normal_weighted.normalized();
        props.centroid = centroid / props.area; // 加权平均
        return props;
    }
    
    void integrateMonomials(const Eigen::MatrixXd& all_nodes, Eigen::MatrixXd& integrals) const override {
        // 后续实现
    }
};