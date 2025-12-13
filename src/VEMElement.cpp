#include "VEMElement.hpp"
#include <iostream>

Eigen::Matrix6d Material::getC() const {
    double l = lambda();
    double m = mu();
    
    Eigen::Matrix6d C = Eigen::Matrix6d::Zero();
    // 对角线
    C(0,0) = C(1,1) = C(2,2) = l + 2*m;
    C(3,3) = C(4,4) = C(5,5) = m;
    
    // 非对角线 (正应力耦合)
    C(0,1) = C(1,0) = C(0,2) = C(2,0) = C(1,2) = C(2,1) = l;
    
    return C;
}

VEMElement::VEMElement(const VEMMesh& mesh, const PolyhedronElement& elem)
    : mesh_(mesh), elem_(elem) {
    // 初始化缩放参数
    // scaling_center_ 取单元质心
    scaling_center_ = elem.centroid;
    
    // scaling_h_ 取单元体积的立方根，或者直径
    // 简单起见，使用 volume^(1/3)
    scaling_h_ = std::pow(elem.volume, 1.0/3.0);
}

Eigen::MatrixXd VEMElement::computeG(const Material& mat) const {
    int n_monos = getNumMonomials();
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n_monos, n_monos);
    Eigen::Matrix6d C = mat.getC();
    double volume = elem_.volume;

    // 定义 helper lambda：计算第 i 个单项式基函数的常数应变 (Voigt notation 6x1)
    // Basis order (total 12):
    // 0-2:  [1,0,0], [0,1,0], [0,0,1]
    // 3-5:  [x~,0,0], [0,x~,0], [0,0,x~]  (x~ = (x-xc)/h)
    // 6-8:  [y~,0,0], [0,y~,0], [0,0,y~]
    // 9-11: [z~,0,0], [0,z~,0], [0,0,z~]
    auto getStrain = [&](int i) -> Eigen::Vector6d {
        Eigen::Vector6d eps = Eigen::Vector6d::Zero();
        double inv_h = 1.0 / scaling_h_;

        // 常数项 (0-2) 的梯度为0，应变也为0
        if (i < 3) return eps;

        // 线性项
        int linear_idx = i - 3; // 0..8
        int dir = linear_idx % 3; // 0=x, 1=y, 2=z (位移分量方向)
        int var = linear_idx / 3; // 0=x, 1=y, 2=z (变量依赖)

        // u = [.., .., ..]
        // du_i / dx_j 
        // strain_xx = du/dx, strain_yy = dv/dy, ...
        // strain_xy = 0.5 * (du/dy + dv/dx)
        
        // 我们的基函数是单分量非零。例如 [x~, 0, 0] -> u=(x-xc)/h, v=0, w=0
        // 只有 du/dx = 1/h 非零。
        
        // 填充 3x3 梯度矩阵 grad(u)
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();
        grad(dir, var) = inv_h; // du_dir / d_var

        // 转换为 Voigt 应变
        // eps = [xx, yy, zz, 2yz, 2xz, 2xy]
        eps(0) = grad(0,0); // xx
        eps(1) = grad(1,1); // yy
        eps(2) = grad(2,2); // zz
        eps(3) = grad(1,2) + grad(2,1); // 2yz
        eps(4) = grad(0,2) + grad(2,0); // 2xz
        eps(5) = grad(0,1) + grad(1,0); // 2xy
        
        return eps;
    };

    // 计算 G_ab = Vol * eps_a^T * C * eps_b
    // 注意：对于 k=1，P1 基函数的应变是常数。所以积分直接提出来变成体积。
    for (int i = 0; i < n_monos; ++i) {
        Eigen::Vector6d eps_i = getStrain(i);
        // 如果应变为0（刚体平移），G 对应的行/列全为 0
        if (eps_i.isZero()) continue;

        for (int j = i; j < n_monos; ++j) {
            Eigen::Vector6d eps_j = getStrain(j);
            if (eps_j.isZero()) continue;

            double val = volume * eps_i.dot(C * eps_j);
            G(i, j) = val;
            G(j, i) = val;
        }
    }

    return G;
}