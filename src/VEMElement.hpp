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
};