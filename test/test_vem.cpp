#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "VEMMesh.hpp"
#include "VEMElement.hpp"

// 辅助打印函数
void checkCondition(bool condition, const std::string& msg) {
    if (condition) {
        std::cout << "[PASSED] " << msg << std::endl;
    } else {
        std::cerr << "[FAILED] " << msg << std::endl;
        exit(1);
    }
}

int main() {
    // 1. 准备网格：标准单位立方体
    VEMMesh mesh;
    Eigen::MatrixXd V(3, 8);
    V << 0, 1, 1, 0, 0, 1, 1, 0,
         0, 0, 1, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 1, 1, 1;
    mesh.setNodes(V);

    // 定义 6 个面
    std::vector<std::vector<int>> faces_nodes = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, // Bottom, Top
        {0, 1, 5, 4}, {1, 2, 6, 5}, // Front, Right
        {2, 3, 7, 6}, {3, 0, 4, 7}  // Back, Left
    };

    std::vector<int> face_ids;
    for (const auto& f : faces_nodes) {
        face_ids.push_back(mesh.addPolygonFace(f));
    }
    mesh.addElement(face_ids);

    // 2. 准备 VEM 计算
    const auto& elem = mesh.getElement(0);
    VEMElement vem_elem(mesh, elem);
    
    Material mat;
    mat.E = 1000.0;
    mat.nu = 0.3;

    // 3. 计算 G 矩阵
    Eigen::MatrixXd G = vem_elem.computeG(mat);
    
    std::cout << "Computed G Matrix size: " << G.rows() << "x" << G.cols() << std::endl;

    // ==========================================
    // 验证 1: 刚体平移 (Translation)
    // 前 3 个基函数是 [1,0,0], [0,1,0], [0,0,1]
    // 它们的应变应为 0，导致 G 的前 3 行/列全为 0
    // ==========================================
    bool trans_zero = true;
    if (G.block(0, 0, 3, 12).norm() > 1e-10) trans_zero = false;
    if (G.block(0, 0, 12, 3).norm() > 1e-10) trans_zero = false;
    checkCondition(trans_zero, "Rigid body translations produce zero energy.");

    // ==========================================
    // 验证 2: 谱分析 (Rank Check)
    // 应该有 6 个特征值近似为 0 (3平移 + 3旋转)
    // ==========================================
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(G);
    Eigen::VectorXd eigenvalues = es.eigenvalues();
    
    int zero_eigenvalues = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-8) { // 稍微放宽一点容差处理浮点误差
            zero_eigenvalues++;
        }
    }
    
    std::cout << "Eigenvalues: " << eigenvalues.transpose() << std::endl;
    checkCondition(zero_eigenvalues == 6, "Matrix G has exactly 6 zero eigenvalues (Rigid Body Modes).");

    // ==========================================
    // 验证 3: 刚体旋转 (Rotation) 检查
    // 构造一个绕 Z 轴旋转的位移场: u = -y, v = x, w = 0
    // 对应基函数系数构造：
    // u ~ -y: dir=0(x), var=1(y) -> linear_idx = 1*3 + 0 = 3 -> index = 3+3 = 6. coeff = -1
    // v ~ x:  dir=1(y), var=0(x) -> linear_idx = 0*3 + 1 = 1 -> index = 1+3 = 4. coeff = 1
    // ==========================================
    Eigen::VectorXd vec_rot_z = Eigen::VectorXd::Zero(12);
    // 注意：需要除以 scaling_h_ 因为基函数是 (x-xc)/h
    // 但为了验证零能，比例系数不影响结果是否为 0
    vec_rot_z(6) = -1.0; 
    vec_rot_z(4) = 1.0;

    double energy_rot_z = vec_rot_z.transpose() * G * vec_rot_z;
    std::cout << "Energy of Z-rotation mode: " << energy_rot_z << std::endl;
    checkCondition(std::abs(energy_rot_z) < 1e-10, "Z-axis rotation produces zero energy.");

    return 0;
}