#include <iostream>
#include <cassert>
#include "VEMMesh.hpp"

// 简单的辅助函数，用于浮点数比较
bool isApprox(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

bool isApproxVec(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double tol = 1e-10) {
    return (a - b).norm() < tol;
}

int main() {
    VEMMesh mesh;

    // 1. 创建节点 (单位立方体)
    Eigen::MatrixXd V(3, 8);
    V << 0, 1, 1, 0, 0, 1, 1, 0,
         0, 0, 1, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 1, 1, 1;
    mesh.setNodes(V);

    // 2. 创建一个底部面 (z=0, 0->3->2->1) 也就是 (0,0)->(0,1)->(1,1)->(1,0)
    // 这是一个正方形，面积应为 1.0，中心在 (0.5, 0.5, 0.0)
    std::vector<int> f1 = {0, 3, 2, 1}; 
    int face_id = mesh.addPolygonFace(f1);

    // 3. 验证单项式积分 (Moments)
    Eigen::Vector4d moments = mesh.getFace(face_id)->integrateMonomials(mesh.getNodes());
    
    std::cout << "Computed Moments: " << moments.transpose() << std::endl;

    // 预期结果
    // m0 (面积) = 1.0
    // m1_x = x_c * Area = 0.5 * 1.0 = 0.5
    // m1_y = y_c * Area = 0.5 * 1.0 = 0.5
    // m1_z = z_c * Area = 0.0 * 1.0 = 0.0
    Eigen::Vector4d expected;
    expected << 1.0, 0.5, 0.5, 0.0;

    if (isApproxVec(moments, expected)) {
        std::cout << "[PASSED] Face moments integration test passed!" << std::endl;
    } else {
        std::cerr << "[FAILED] Face moments mismatch!" << std::endl;
        std::cerr << "Expected: " << expected.transpose() << std::endl;
        return 1;
    }

    return 0;
}