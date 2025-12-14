#pragma once
#include <cmath>
#include <functional>

namespace vem::utils {

/**
 * @brief 欧拉-伯努利悬臂梁理论解
 * * 假设梁固定在 x_fixed，向 x 正方向延伸。
 * 重力/载荷 q 沿 Z 轴负方向。
 */
class CantileverBeamSolution {
public:
    double E; // 杨氏模量
    double I; // 截面惯性矩 (w*h^3/12)
    double q; // 线载荷 (rho * Area * g)
    double L; // 梁长
    double x_fixed; // 固定端 X 坐标

    CantileverBeamSolution(double E, double I, double q, double L, double x_fixed = -2.0)
        : E(E)
        , I(I)
        , q(q)
        , L(L)
        , x_fixed(x_fixed)
    {
    }

    // 获取 (x,y,z) 处的理论 Z 向位移
    double getDisplacementZ(double x) const
    {
        double x_loc = x - x_fixed; // 局部坐标 [0, L]
        if (x_loc < 0)
            x_loc = 0;
        if (x_loc > L)
            x_loc = L;

        // w(x) = - (q x^2) / (24 EI) * (x^2 - 4Lx + 6L^2)
        // 负号表示沿 Z 轴向下
        double w = -(q * x_loc * x_loc) / (24.0 * E * I) * (x_loc * x_loc - 4.0 * L * x_loc + 6.0 * L * L);
        return w;
    }
};

} // namespace vem::utils
