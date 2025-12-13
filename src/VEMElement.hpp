#pragma once
#include <Eigen/Dense>
#include "VEMMesh.hpp"  

namespace Eigen
{
    using Matrix6d = Eigen::Matrix<double, 6, 6>;
    using Vector6d = Eigen::Matrix<double, 6, 1>;
} // namespace Eigen



// 简单的线性弹性材料参数
struct Material {
    double E;  // Young's modulus
    double nu; // Poisson's ratio

    // 获取 Lame 常数 lambda
    double lambda() const {
        return (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
    }

    // 获取 Lame 常数 mu (Shear modulus)
    double mu() const {
        return E / (2.0 * (1.0 + nu));
    }
    
    // 计算 Voigt 表示下的弹性张量 C (6x6 矩阵)
    // sigma = C * epsilon
    // order: xx, yy, zz, yz, xz, xy
    Eigen::Matrix6d getC() const;
};

class VEMElement {
private:
    const VEMMesh& mesh_;
    const PolyhedronElement& elem_;
    
    // 缩放参数，用于提高条件数
    double scaling_h_; 
    Eigen::Vector3d scaling_center_;

public:
    VEMElement(const VEMMesh& mesh, const PolyhedronElement& elem);

    // 计算一致性矩阵 G (Matrix G)
    // G_ab = int_E strain(m_a) : C : strain(m_b) dV
    // 对于 k=1，应变是常数，积分简化为 Volume * ...
    Eigen::MatrixXd computeG(const Material& mat) const;
    
    // 辅助：获取线性弹性的一阶单项式基函数的数量 (3D向量场)
    // k=1 时，基函数包含：3个刚体平移 + 3个线性变换 * 3个方向 = 12 个基函数？
    // Gain 2014 使用的基是多项式空间 P1 的基。
    // 基函数顺序建议：
    // [1, 0, 0], [0, 1, 0], [0, 0, 1], (平移)
    // [(x-xc)/h, 0, 0], [0, (x-xc)/h, 0], ... 
    static int getNumMonomials() { return 12; } // 3 (constant) + 9 (linear)

    // 计算右端项矩阵 B (Matrix B)
    // B_aI = int_dE (C : eps(m_a) . n) . phi_I dS
    // Returns matrix of size (n_monos x n_dofs)
    Eigen::MatrixXd computeB(const Material& mat) const;

    // 辅助：获取自由度总数 (3 * N_vertices)
    int getNumDofs() const { return 3 * elem_.face_indices.size(); /* 注意：这不对，elem 只有 face_indices */ }
    // 修正：我们需要从 mesh 获取单元的顶点列表。
    // 建议在 PolyhedronElement 中存储顶点列表，或者在此处实时收集。
    // 为了简单，我们可以在 computeB 内部收集顶点。

    // 新增：计算 D 矩阵 (Matrix D)
    // D_i_alpha = value of monomial alpha at DOF i
    // Size: (n_dofs x n_monos)
    Eigen::MatrixXd computeD() const;

    // 新增：计算单元刚度矩阵 K = Kc + Ks
    Eigen::MatrixXd computeStiffness(const Material& mat) const;
};