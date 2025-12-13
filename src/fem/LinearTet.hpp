#pragma once
#include <Eigen/Dense>
#include "../VEMElement.hpp" // 复用 Material 定义

namespace vem::fem {

class LinearTet {
public:
    // 计算 4 节点四面体的刚度矩阵 (12x12)
    // nodes: 3x4 矩阵，每一列是一个顶点的坐标
    static Eigen::MatrixXd computeStiffness(const Eigen::MatrixXd& nodes, const Material& mat);

    // 辅助：计算 B 矩阵 (6x12) 和 体积
    // 返回 pair<B matrix, Volume>
    static std::pair<Eigen::MatrixXd, double> computeB_and_Volume(const Eigen::MatrixXd& nodes);
};

} // namespace vem::fem