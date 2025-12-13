#include "LinearTet.hpp"
#include <iostream>

namespace vem::fem {

std::pair<Eigen::MatrixXd, double> LinearTet::computeB_and_Volume(const Eigen::MatrixXd& nodes) {
    // 1. 构建 P 矩阵 (4x4)
    Eigen::Matrix4d P;
    P.row(0).setOnes();
    P.block(1, 0, 3, 4) = nodes; // Row 1..3 = x,y,z

    double detP = P.determinant();
    double volume = std::abs(detP) / 6.0;

    Eigen::Matrix4d P_inv = P.inverse();
    
    // 2. 提取梯度 dN
    // P_inv 的每一行 i 对应节点 i 的系数 [a, b, c, d]
    // 梯度 grad(Ni) = [b, c, d]^T，位于 P_inv 的第 1,2,3 列
    
    // 我们需要 dN 矩阵 (3x4)，其中第 i 列是节点 i 的梯度
    // P_inv.block(0, 1, 4, 3) 取出的是 [4 nodes x 3 coeffs]
    // 所以需要转置
    Eigen::MatrixXd dN = P_inv.block(0, 1, 4, 3).transpose(); 

    // 3. 构建 B 矩阵 (6x12) - Voigt Notation
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 12);
    
    for (int i = 0; i < 4; ++i) {
        double b = dN(0, i); // dNi/dx
        double c = dN(1, i); // dNi/dy
        double d = dN(2, i); // dNi/dz
        
        int col = 3 * i;
        
        B(0, col) = b;   // xx
        B(1, col+1) = c; // yy
        B(2, col+2) = d; // zz
        
        B(3, col+1) = d; B(3, col+2) = c; // yz
        B(4, col) = d;   B(4, col+2) = b; // xz
        B(5, col) = c;   B(5, col+1) = b; // xy
    }

    return {B, volume};
}

Eigen::MatrixXd LinearTet::computeStiffness(const Eigen::MatrixXd& nodes, const Material& mat) {
    auto [B, volume] = computeB_and_Volume(nodes);
    Eigen::Matrix6d C = mat.getC();
    return volume * B.transpose() * C * B;
}

} // namespace vem::fem