#include "Geometry.hpp"



    GeometricProps TriFace::computeProps(const Eigen::MatrixXd& all_nodes) const {
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
    void TriFace::integrateMonomials(const Eigen::MatrixXd& all_nodes, Eigen::MatrixXd& integrals) const{

    }

// 在 Geometry.cpp 中补充 PolygonFace 的实现

GeometricProps PolygonFace::computeProps(const Eigen::MatrixXd& all_nodes) const {
    GeometricProps props;
    props.area = 0.0;
    props.centroid.setZero();
    props.normal.setZero();

    // 1. 计算均值中心作为临时参考点 (避免凹多边形甚至某些复杂情况的误差，但对凸多边形最稳)
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int idx : node_indices) center += all_nodes.col(idx);
    center /= node_indices.size();

    // 2. 扇形剖分 (Fan triangulation) 累加面积和一阶矩
    // 将多边形分解为 [center, node_i, node_{i+1}] 的三角形集合
    for (size_t i = 0; i < node_indices.size(); ++i) {
        Eigen::Vector3d p1 = all_nodes.col(node_indices[i]);
        Eigen::Vector3d p2 = all_nodes.col(node_indices[(i + 1) % node_indices.size()]); // 循环连接

        Eigen::Vector3d edge1 = p1 - center;
        Eigen::Vector3d edge2 = p2 - center;
        Eigen::Vector3d cross = edge1.cross(edge2); // 局部三角形法向 * 2倍面积

        double tri_area = 0.5 * cross.norm();
        Eigen::Vector3d tri_centroid = (center + p1 + p2) / 3.0;

        // 累加
        props.area += tri_area;
        props.centroid += tri_centroid * tri_area;
        props.normal += cross; // 注意：这里简单累加可能需要根据多边形共面性调整，通常取平均或第一个有效法向
    }

    if (props.area > 1e-12) {
        props.centroid /= props.area;
        props.normal.normalize();
    }
    return props;
}