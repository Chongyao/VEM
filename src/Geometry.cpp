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