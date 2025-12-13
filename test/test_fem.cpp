#include <iostream>
#include "VEMMesh.hpp"
#include "fem/FEMSolver.hpp"

// 线性解析解
Eigen::Vector3d exactSolution(const Eigen::Vector3d& p) {
    return Eigen::Vector3d(0.1 * p.x(), 0.0, 0.0);
}

int main() {
    // 1. 创建一个单四面体网格
    VEMMesh mesh;
    Eigen::MatrixXd V(3, 4);
    V << 0, 1, 0, 0,
         0, 0, 1, 0,
         0, 0, 0, 1;
    mesh.setNodes(V);
    
    // Tet faces: (0,2,1), (0,1,3), (0,3,2), (1,2,3)
    // 只要是闭合的即可
    std::vector<int> f1 = {0, 2, 1};
    std::vector<int> f2 = {0, 1, 3};
    std::vector<int> f3 = {0, 3, 2};
    std::vector<int> f4 = {1, 2, 3};
    
    std::vector<int> face_ids;
    face_ids.push_back(mesh.addTriFace(f1));
    face_ids.push_back(mesh.addTriFace(f2));
    face_ids.push_back(mesh.addTriFace(f3));
    face_ids.push_back(mesh.addTriFace(f4));
    
    mesh.addElement(face_ids);
    
    // 2. FEM 求解
    vem::fem::FEMSolver solver(mesh);
    Material mat; mat.E = 1000; mat.nu = 0.3;
    solver.assemble(mat);
    
    // 3. Patch Test (固定 3 个点，解第 4 个点)
    std::map<int, double> bc;
    for(int i=0; i<3; ++i) { // Fix nodes 0,1,2
        Eigen::Vector3d u = exactSolution(mesh.getNodes().col(i));
        bc[3*i+0]=u(0); bc[3*i+1]=u(1); bc[3*i+2]=u(2);
    }
    // Node 3 is free
    
    Eigen::VectorXd f = Eigen::VectorXd::Zero(12);
    Eigen::VectorXd u_sol = solver.solve(bc, f);
    
    // 4. 验证
    Eigen::Vector3d u3 = u_sol.segment<3>(3*3);
    Eigen::Vector3d u3_ex = exactSolution(mesh.getNodes().col(3));
    
    std::cout << "FEM Node 3: " << u3.transpose() << std::endl;
    std::cout << "Exact     : " << u3_ex.transpose() << std::endl;
    
    if ((u3 - u3_ex).norm() < 1e-10) {
        std::cout << "[PASSED] FEM Linear Tet Test" << std::endl;
        return 0;
    } else {
        std::cout << "[FAILED] FEM Linear Tet Test" << std::endl;
        return 1;
    }
}