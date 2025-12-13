#include "LinearTet.hpp"
#include <iostream>

namespace vem::fem {

std::pair<Eigen::MatrixXd, double> LinearTet::computeB_and_Volume(const Eigen::MatrixXd& nodes) {
    // 1. 构建形状函数矩阵 P
    // [1  1  1  1 ]
    // [x1 x2 x3 x4]
    // [y1 y2 y3 y4]
    // [z1 z2 z3 z4]
    Eigen::Matrix4d P;
    P.row(0).setOnes();
    P.block(1, 0, 3, 4) = nodes;

    // 2. 求逆获得形状函数系数
    // P * Coeffs = I  =>  Coeffs = P^-1
    // N_i = a_i + b_i*x + c_i*y + d_i*z
    // Coeffs 矩阵的每一列对应 [a, b, c, d]_i ? 
    // 不，通常 N = [1 x y z] * C. 
    // [1 x_j y_j z_j] * C_col_i = delta_ji
    // 所以 P.transpose() * C = I  => C = P^-T
    
    // 简单方法：计算梯度 Grad N
    // Grad N 是常数矩阵 (3x4)
    // Grad N = [dN1/dx ...; dN1/dy ...; dN1/dz ...]
    
    double detP = P.determinant();
    double volume = std::abs(detP) / 6.0;

    // P_inv 的行包含了梯度信息
    // [ a1 a2 a3 a4 ]
    // [ b1 b2 b3 b4 ] -> dN/dx
    // [ c1 c2 c3 c4 ] -> dN/dy
    // [ d1 d2 d3 d4 ] -> dN/dz
    Eigen::Matrix4d P_inv = P.inverse();
    
    // 提取梯度 (后三行)
    // 每一列对应一个节点：Col 0 -> Node 0
    Eigen::MatrixXd dN = P_inv.block(1, 0, 3, 4); // 3x4 Matrix: [dN_i / dx_j]

    // 3. 构建 B 矩阵 (6x12) - Voigt Notation
    // Strain = B * u
    // u = [u1 v1 w1 u2 v2 w2 ...]^T
    // eps = [xx yy zz 2yz 2xz 2xy]^T (注意 Eigen 习惯通常工程剪切应变包含 2)
    // Material::getC 假设的是 Voigt: xx, yy, zz, yz, xz, xy.
    // 所以我们的 B 矩阵要匹配 Material 的定义。
    // 注意：Material::getC 里通常是 lambda*(trace) + 2mu*eps. 
    // 如果 Material 用标准工程应变 (gamma = 2*eps_tensor)，则 B 矩阵对应行应该是 dNi/dy + dNj/dx.
    // 让我们假设 computeG 中的逻辑：
    // sigma = C * eps. 
    // 这里的 eps 是 tensor component 还是 engineering? 
    // 通常 linear elasticity FEM assembly 使用 engineering strains [xx, yy, zz, g_yz, g_xz, g_xy]
    // 你的 computeG 里： eps(3) = grad(1,2) + grad(2,1); // 2yz = gamma_yz
    // 所以是的，我们需要 Engineering Strain。

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 12);
    
    for (int i = 0; i < 4; ++i) {
        double b = dN(0, i); // dNi/dx
        double c = dN(1, i); // dNi/dy
        double d = dN(2, i); // dNi/dz
        
        int col = 3 * i;
        
        // xx
        B(0, col) = b;
        // yy
        B(1, col+1) = c;
        // zz
        B(2, col+2) = d;
        // yz (gamma) = dv/dz + dw/dy
        B(3, col+1) = d; B(3, col+2) = c;
        // xz (gamma) = du/dz + dw/dx
        B(4, col) = d;   B(4, col+2) = b;
        // xy (gamma) = du/dy + dv/dx
        B(5, col) = c;   B(5, col+1) = b;
    }

    return {B, volume};
}

Eigen::MatrixXd LinearTet::computeStiffness(const Eigen::MatrixXd& nodes, const Material& mat) {
    auto [B, volume] = computeB_and_Volume(nodes);
    Eigen::Matrix6d C = mat.getC();
    
    // K = integral B^T C B dV = Vol * B^T * C * B
    return volume * B.transpose() * C * B;
}

} // namespace vem::fem