#include "Geometry.hpp"

// --- TriFace 实现 ---

GeometricProps TriFace::computeProps(const Eigen::MatrixXd& all_nodes) const {
    Eigen::Vector3d p0 = all_nodes.col(node_indices[0]);
    Eigen::Vector3d p1 = all_nodes.col(node_indices[1]);
    Eigen::Vector3d p2 = all_nodes.col(node_indices[2]);

    Eigen::Vector3d cross = (p1 - p0).cross(p2 - p0);
    
    GeometricProps props;
    props.area = 0.5 * cross.norm();
    props.normal = cross.normalized();
    props.centroid = (p0 + p1 + p2) / 3.0;
    return props;
}

// --- TriFace Clone ---
std::unique_ptr<VEMFace> TriFace::clone() const {
    // 调用默认拷贝构造函数
    return std::make_unique<TriFace>(*this);
}

Eigen::Vector4d TriFace::integrateMonomials(const Eigen::MatrixXd& all_nodes) const {
    // 复用 computeProps 的逻辑，减少代码重复
    GeometricProps props = computeProps(all_nodes);
    
    Eigen::Vector4d moments;
    moments(0) = props.area; // 0阶矩: 面积
    // 1阶矩: 质心 * 面积 (即 int x dA, int y dA, int z dA)
    moments.segment<3>(1) = props.centroid * props.area;
    
    return moments;
}

// --- PolygonFace 实现 ---

GeometricProps PolygonFace::computeProps(const Eigen::MatrixXd& all_nodes) const {
    GeometricProps props;
    props.area = 0.0;
    props.centroid.setZero();
    props.normal.setZero();

    // 计算均值中心
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int idx : node_indices) center += all_nodes.col(idx);
    center /= node_indices.size();

    // 扇形剖分累加
    for (size_t i = 0; i < node_indices.size(); ++i) {
        Eigen::Vector3d p1 = all_nodes.col(node_indices[i]);
        Eigen::Vector3d p2 = all_nodes.col(node_indices[(i + 1) % node_indices.size()]);

        Eigen::Vector3d cross = (p1 - center).cross(p2 - center);
        double tri_area = 0.5 * cross.norm();
        Eigen::Vector3d tri_centroid = (center + p1 + p2) / 3.0;

        props.area += tri_area;
        props.centroid += tri_centroid * tri_area;
        props.normal += cross; 
    }

    if (props.area > 1e-12) {
        props.centroid /= props.area;
        props.normal.normalize();
    }
    return props;
}

Eigen::Vector4d PolygonFace::integrateMonomials(const Eigen::MatrixXd& all_nodes) const {
    // 对于一阶 VEM，面的单项式积分本质上就是计算面积和质心矩
    // 所以这里的逻辑和 computeProps 几乎一样，只是返回形式不同。
    // 在高阶 VEM (k>=2) 时，这里需要计算二阶矩等，代码会变得不同。
    
    GeometricProps props = computeProps(all_nodes);
    
    Eigen::Vector4d moments;
    moments(0) = props.area;
    moments.segment<3>(1) = props.centroid * props.area;
    
    return moments;
}

// --- TriFace 新增实现 ---
std::map<int, double> TriFace::computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const {
    // 三角形的形状函数积分非常简单：Area / 3
    GeometricProps props = computeProps(all_nodes);
    double weight = props.area / 3.0;
    
    std::map<int, double> weights;
    for (int idx : node_indices) {
        weights[idx] = weight;
    }
    return weights;
}

// --- PolygonFace 新增实现 ---
std::map<int, double> PolygonFace::computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const {
    std::map<int, double> weights;
    // 初始化权重为0
    for (int idx : node_indices) weights[idx] = 0.0;

    // 1. 计算中心
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int idx : node_indices) center += all_nodes.col(idx);
    center /= node_indices.size();

    double total_area = 0.0;
    int N = node_indices.size();

    // 2. 扇形剖分积分
    // 假设 VEM 面内形状函数使得：中心点的值 v(C) = 1/N * sum(v_i)
    // 且在每个子三角形上是线性的。
    // 积分 int_f phi_i = sum_k int_{Tk} phi_i
    
    for (size_t k = 0; k < N; ++k) {
        int idx_a = node_indices[k];
        int idx_b = node_indices[(k + 1) % N];
        
        Eigen::Vector3d pa = all_nodes.col(idx_a);
        Eigen::Vector3d pb = all_nodes.col(idx_b);

        // 子三角形 T_k = (Center, A, B)
        double tri_area = 0.5 * (pa - center).cross(pb - center).norm();
        total_area += tri_area;

        // 在三角形 (C, A, B) 上的积分贡献：
        // 节点 A 的贡献: tri_area / 3
        // 节点 B 的贡献: tri_area / 3
        // 节点 Center 的贡献: tri_area / 3
        
        weights[idx_a] += tri_area / 3.0;
        weights[idx_b] += tri_area / 3.0;
        
        // Center 的贡献均分给所有顶点 (因为 v(C) = mean(v_i))
        double center_contrib = (tri_area / 3.0) / N;
        for (int idx : node_indices) {
            weights[idx] += center_contrib;
        }
    }
    
    return weights;
}

std::unique_ptr<VEMFace> PolygonFace::clone() const {
    return std::make_unique<PolygonFace>(*this);
}