#include <iostream>
#include <iomanip>
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"
#include "fem/FEMSolver.hpp"

// 辅助：计算矩阵的相对误差
double computeMatrixError(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    double normA = A.norm();
    if (normA < 1e-12) return (A - B).norm(); // 避免除零
    return (A - B).norm() / normA;
}

int main() {
    std::cout << "=== VEM vs FEM Equivalence Test (on Tetrahedron) ===" << std::endl;

    // 1. 创建单四面体网格
    // 为了更有代表性，我们用一个稍微不规则的四面体
    VEMMesh mesh;
    Eigen::MatrixXd V(3, 4);
    V << 0.0, 1.2, 0.1, 0.2,  // x
         0.0, 0.1, 1.1, 0.3,  // y
         0.0, 0.1, 0.2, 1.3;  // z
    mesh.setNodes(V);

    std::vector<int> face_ids;
    face_ids.push_back(mesh.addTriFace({0, 2, 1}));
    face_ids.push_back(mesh.addTriFace({0, 1, 3}));
    face_ids.push_back(mesh.addTriFace({0, 3, 2}));
    face_ids.push_back(mesh.addTriFace({1, 2, 3}));
    mesh.addElement(face_ids);

    Material mat; 
    mat.E = 10000.0; 
    mat.nu = 0.25;

    // 2. 计算 FEM 刚度矩阵
    vem::fem::FEMSolver fem_solver(mesh);
    fem_solver.assemble(mat);
    Eigen::MatrixXd K_fem = Eigen::MatrixXd(fem_solver.getK()); // Convert sparse to dense

    // 3. 计算 VEM 刚度矩阵
    VEMSolver vem_solver(mesh);
    vem_solver.assemble(mat);
    Eigen::MatrixXd K_vem = Eigen::MatrixXd(vem_solver.getK());

    // 4. 对比
    // 由于浮点数运算顺序不同（VEM 涉及大量投影、求逆），
    // 我们期望误差在机器精度允许范围内（例如 1e-10 ~ 1e-13），而不是完全的二进制 0。
    
    double error = computeMatrixError(K_vem, K_fem);
    
    std::cout << "K_fem norm: " << K_fem.norm() << std::endl;
    std::cout << "K_vem norm: " << K_vem.norm() << std::endl;
    std::cout << "Relative Error ||K_vem - K_fem|| / ||K_fem||: " << error << std::endl;

    // 打印前几行对比 (Debugging)
    std::cout << "\nComparison (Top-Left 6x6 block):" << std::endl;
    std::cout << "FEM:\n" << K_fem.block(0,0,6,6) << "\n" << std::endl;
    std::cout << "VEM:\n" << K_vem.block(0,0,6,6) << "\n" << std::endl;
    std::cout << "Diff:\n" << (K_vem - K_fem).block(0,0,6,6) << std::endl;

    if (error < 1e-10) {
        std::cout << "\n[PASSED] VEM is equivalent to FEM on Tetrahedra!" << std::endl;
        return 0;
    } else {
        std::cout << "\n[FAILED] VEM differs from FEM too much." << std::endl;
        return 1;
    }
}