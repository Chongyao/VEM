#pragma once
#include "VEMMesh.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>

namespace vem::utils {

/**
 * @brief 计算节点位移的相对 L2 误差
 * * Error = || u_num - u_exact ||_2 / || u_exact ||_2
 * * @param mesh 网格对象
 * @param u_vem VEM 计算出的位移向量 (大小 3*N)
 * @param exact_func 精确解函数，输入 (x,y,z)，返回 vec3 位移
 * @return double 相对误差
 */
inline double computeDiscreteL2Error(const VEMMesh& mesh,
    const Eigen::VectorXd& u_vem,
    std::function<Eigen::Vector3d(double, double, double)> exact_func)
{
    double error_sq_sum = 0.0;
    double ref_sq_sum = 0.0;

    const auto& nodes = mesh.getNodes();
    int n_nodes = mesh.getNumNodes();

    for (int i = 0; i < n_nodes; ++i) {
        double x = nodes(0, i);
        double y = nodes(1, i);
        double z = nodes(2, i);

        // 数值解
        Eigen::Vector3d u_num = u_vem.segment<3>(3 * i);

        // 精确解
        Eigen::Vector3d u_ref = exact_func(x, y, z);

        // 累加
        error_sq_sum += (u_num - u_ref).squaredNorm();
        ref_sq_sum += u_ref.squaredNorm();
    }

    if (ref_sq_sum < 1e-15)
        return 0.0; // 避免除零
    return std::sqrt(error_sq_sum / ref_sq_sum);
}

} // namespace vem::utils
