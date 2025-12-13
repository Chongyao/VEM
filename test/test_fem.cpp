#include <iostream>
#include "VEMMesh.hpp"
#include "fem/FEMSolver.hpp"
#include "fem/LinearTet.hpp" // 为了计算理论力

// 线性解析解 u(x) = (0.1x, 0, 0)
Eigen::Vector3d exactSolution(const Eigen::Vector3d& p) {
    return Eigen::Vector3d(0.1 * p.x(), 0.0, 0.0);
}

int main() {
    // 1. 创建单四面体
    VEMMesh mesh;
    Eigen::MatrixXd V(3, 4);
    V << 0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;
    mesh.setNodes(V);
    
    std::vector<int> face_ids;
    face_ids.push_back(mesh.addTriFace({0, 2, 1}));
    face_ids.push_back(mesh.addTriFace({0, 1, 3}));
    face_ids.push_back(mesh.addTriFace({0, 3, 2}));
    face_ids.push_back(mesh.addTriFace({1, 2, 3}));
    mesh.addElement(face_ids);
    
    // 2. 准备求解器
    vem::fem::FEMSolver solver(mesh);
    Material mat; mat.E = 1000; mat.nu = 0.3;
    solver.assemble(mat);
    
    // 3. 边界条件
    // Dirichlet: 固定节点 0, 1, 2 (底面)
    std::map<int, double> bc;
    for(int i=0; i<3; ++i) { 
        Eigen::Vector3d u = exactSolution(mesh.getNodes().col(i));
        bc[3*i+0]=u(0); bc[3*i+1]=u(1); bc[3*i+2]=u(2);
    }
    
    // 4. 计算并施加理论外力 (Neumann BC) 到自由节点 3
    // 只有这样，物理上才能维持非物理的单轴应变状态 (无泊松收缩)
    // f_ext = Integral(B^T * sigma) dV
    // sigma = C * epsilon_exact
    
    Eigen::Vector6d eps_exact = Eigen::Vector6d::Zero();
    eps_exact(0) = 0.1; // xx strain
    Eigen::Vector6d sigma_exact = mat.getC() * eps_exact;
    
    auto [B, Vol] = vem::fem::LinearTet::computeB_and_Volume(V);
    Eigen::VectorXd f_elem = Vol * B.transpose() * sigma_exact;
    
    Eigen::VectorXd f_global = Eigen::VectorXd::Zero(12);
    // 将单元力加到全局力 (这里单元节点序就是全局序 0,1,2,3)
    f_global = f_elem;
    
    // 5. 求解
    Eigen::VectorXd u_sol = solver.solve(bc, f_global);
    
    // 6. 验证
    Eigen::Vector3d u3 = u_sol.segment<3>(3*3);
    Eigen::Vector3d u3_ex = exactSolution(mesh.getNodes().col(3));
    
    std::cout << "FEM Node 3: " << u3.transpose() << std::endl;
    std::cout << "Exact     : " << u3_ex.transpose() << std::endl;
    
    if ((u3 - u3_ex).norm() < 1e-10) {
        std::cout << "[PASSED] FEM Linear Tet Test (with Consistent Forces)" << std::endl;
        return 0;
    } else {
        std::cout << "[FAILED] FEM Linear Tet Test" << std::endl;
        return 1;
    }
}