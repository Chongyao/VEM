#pragma once
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <iostream>

namespace vem {

class MathUtils {
public:
    /**
     * @brief 计算刚度矩阵的稳定性比率 sigma = lambda_min / lambda_max
     * 忽略接近 0 的刚体模态
     * @param K 单元刚度矩阵
     * @return double 稳定性比率 (0.0 ~ 1.0)
     */
    static double computeStabilityRatio(const Eigen::MatrixXd& K)
    {
        // 1. 计算特征值 (SelfAdjoint 用于对称矩阵，速度更快)
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
        const auto& eigenvalues = es.eigenvalues();

        double min_lambda = 1e30;
        double max_lambda = -1e30;

        // 2. 确定最大特征值
        double max_val = eigenvalues.maxCoeff();
        // 如果矩阵全零或极小，直接返回 0 (极度病态)
        if (max_val < 1e-15)
            return 0.0;

        // 3. 设定阈值过滤刚体模态 (Rigid Body Modes)
        // 3D 弹性力学理论上有 6 个零能模式 (3平移+3旋转)
        // 我们设定一个相对阈值，比如最大特征值的 1e-8
        double threshold = max_val * 1e-8;

        for (int i = 0; i < eigenvalues.size(); ++i) {
            double val = eigenvalues(i);
            if (val > max_lambda)
                max_lambda = val;

            // 只统计非零特征值
            if (val > threshold && val < min_lambda)
                min_lambda = val;
        }

        // 如果没找到非零特征值（比如全零矩阵），返回 0
        if (min_lambda == 1e30)
            return 0.0;

        return min_lambda / max_lambda;
    }
};

} // namespace vem
