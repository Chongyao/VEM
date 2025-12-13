#include <iostream>
#include <vector>
#include <map>
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"

// 任意线性位移场
// u = a0 + a1*x + a2*y + a3*z
Eigen::Vector3d exactSolution(const Eigen::Vector3d& p) {
    return Eigen::Vector3d(
        0.1 + 0.2 * p.x() + 0.1 * p.y() + 0.05 * p.z(),
        0.05 - 0.1 * p.x() + 0.2 * p.y() - 0.1 * p.z(),
        -0.1 + 0.05 * p.x() + 0.1 * p.y() + 0.3 * p.z()
    );
}

int main() {
    std::cout << "=== 2x2x2 Internal Node Patch Test ===" << std::endl;

    // 1. 构建 2x2x2 网格 (8个单元, 27个节点)
    // 节点网格: x,y,z in {0, 1, 2}
    VEMMesh mesh;
    int nx=2, ny=2, nz=2;
    int n_nodes = 27; // 3*3*3
    Eigen::MatrixXd V(3, n_nodes);
    
    int internal_node_id = -1;

    for(int k=0; k<=nz; ++k) {
        for(int j=0; j<=ny; ++j) {
            for(int i=0; i<=nx; ++i) {
                int idx = k*(nx+1)*(ny+1) + j*(nx+1) + i;
                V(0, idx) = i * 1.0; // dx=1
                V(1, idx) = j * 1.0;
                V(2, idx) = k * 1.0;
                
                // 找到中心点 (1,1,1)
                if (i==1 && j==1 && k==1) internal_node_id = idx;
            }
        }
    }
    mesh.setNodes(V);
    std::cout << "Internal Node ID: " << internal_node_id << std::endl;

    // 简单的六面体生成器
    auto get_node = [&](int i, int j, int k) {
        return k*(nx+1)*(ny+1) + j*(nx+1) + i;
    };

    auto add_hex = [&](int i, int j, int k) {
        // Hexahedron nodes
        int n0 = get_node(i,j,k);     int n1 = get_node(i+1,j,k);
        int n2 = get_node(i+1,j+1,k); int n3 = get_node(i,j+1,k);
        int n4 = get_node(i,j,k+1);   int n5 = get_node(i+1,j,k+1);
        int n6 = get_node(i+1,j+1,k+1); int n7 = get_node(i,j+1,k+1);
        
        std::vector<int> faces;
        // 顺序很重要 (Outward normal)
        faces.push_back(mesh.addPolygonFace({n0, n3, n2, n1})); // Bottom
        faces.push_back(mesh.addPolygonFace({n4, n5, n6, n7})); // Top
        faces.push_back(mesh.addPolygonFace({n0, n1, n5, n4})); // Front
        faces.push_back(mesh.addPolygonFace({n1, n2, n6, n5})); // Right
        faces.push_back(mesh.addPolygonFace({n2, n3, n7, n6})); // Back
        faces.push_back(mesh.addPolygonFace({n3, n0, n4, n7})); // Left
        
        mesh.addElement(faces);
    };

    for(int k=0; k<nz; ++k)
        for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
                add_hex(i, j, k);

    // 2. 求解设置
    VEMSolver solver(mesh);
    Material mat; mat.E = 1000; mat.nu = 0.3;
    solver.assemble(mat);

    // 3. 边界条件: 固定除中心点以外的所有节点
    std::map<int, double> dirichlet_bc;
    for (int i = 0; i < n_nodes; ++i) {
        if (i == internal_node_id) continue; // 让中心点自由
        
        Eigen::Vector3d u_ex = exactSolution(mesh.getNodes().col(i));
        dirichlet_bc[3*i+0] = u_ex(0);
        dirichlet_bc[3*i+1] = u_ex(1);
        dirichlet_bc[3*i+2] = u_ex(2);
    }
    
    Eigen::VectorXd f_ext = Eigen::VectorXd::Zero(solver.getNumDofs()); // 无外力

    // 4. 求解
    Eigen::VectorXd u_sol = solver.solve(dirichlet_bc, f_ext);
    
    // 5. 验证中心点
    Eigen::Vector3d u_internal = u_sol.segment<3>(3*internal_node_id);
    Eigen::Vector3d u_exact = exactSolution(mesh.getNodes().col(internal_node_id));
    
    std::cout << "Internal Node Pos: " << mesh.getNodes().col(internal_node_id).transpose() << std::endl;
    std::cout << "Computed Disp: " << u_internal.transpose() << std::endl;
    std::cout << "Exact    Disp: " << u_exact.transpose() << std::endl;
    
    double err = (u_internal - u_exact).norm();
    std::cout << "Error: " << err << std::endl;
    
    if (err < 1e-10) {
        std::cout << "[PASSED] Correct Patch Test Passed!" << std::endl;
        return 0;
    } else {
        std::cout << "[FAILED] Patch Test Failed." << std::endl;
        return 1;
    }
}