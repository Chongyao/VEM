#include <iostream>
#include <vector>
#include <cmath>
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"

// 辅助: 线性解析解 u(x,y,z)
// 我们可以设一个比较复杂的线性场，不仅是常数
// u = 0.1 + 0.2x - 0.05y + 0.1z
// v = ...
Eigen::Vector3d exactSolution(const Eigen::Vector3d& p) {
    Eigen::Vector3d u;
    u(0) = 0.1 + 0.2 * p.x() + 0.1 * p.y() + 0.05 * p.z();
    u(1) = 0.05 - 0.1 * p.x() + 0.2 * p.y() - 0.1 * p.z();
    u(2) = -0.1 + 0.0 * p.x() + 0.1 * p.y() + 0.3 * p.z();
    return u;
}

int main() {
    // 1. 创建一个 2x2x2 的简单网格 (8个单元)
    VEMMesh mesh;
    
    int nx = 2, ny = 2, nz = 2;
    int n_nodes = (nx+1)*(ny+1)*(nz+1);
    Eigen::MatrixXd V(3, n_nodes);
    
    // 生成规则网格点
    int idx = 0;
    for(int k=0; k<=nz; ++k)
        for(int j=0; j<=ny; ++j)
            for(int i=0; i<=nx; ++i) {
                V(0, idx) = i; // dx=1
                V(1, idx) = j;
                V(2, idx) = k;
                idx++;
            }
    mesh.setNodes(V);
    
    // 生成六面体单元 (这里简化，假设我们手动或用工具生成了 connectivity)
    // 为了简单起见，我们只用刚才 test_solver 的 2 单元逻辑，或者稍微扩展
    // 让我们复用 test_solver 的 2 单元网格，这足够做 Patch Test 了
    // 重新定义网格: 2个并排立方体
    mesh = VEMMesh(); // Reset
    Eigen::MatrixXd V2(3, 12);
    V2 << 0,0,0,0, 1,1,1,1, 2,2,2,2, // x
          0,1,1,0, 0,1,1,0, 0,1,1,0, // y
          0,0,1,1, 0,0,1,1, 0,0,1,1; // z
    mesh.setNodes(V2);
    
    auto addFace = [&](int n1, int n2, int n3, int n4) {
        return mesh.addPolygonFace({n1, n2, n3, n4});
    };
    int f1_back   = addFace(0, 3, 2, 1);
    int f1_front  = addFace(4, 5, 6, 7);
    int f1_bottom = addFace(0, 1, 5, 4);
    int f1_top    = addFace(3, 7, 6, 2);
    int f1_left   = addFace(0, 4, 7, 3);
    int f1_right  = addFace(1, 2, 6, 5);
    mesh.addElement({f1_back, f1_front, f1_bottom, f1_top, f1_left, f1_right});

    int f2_front  = addFace(8, 9, 10, 11);
    int f2_bottom = addFace(4, 5, 9, 8);
    int f2_top    = addFace(7, 11, 10, 6);
    int f2_left   = addFace(4, 8, 11, 7);
    int f2_right  = addFace(5, 6, 10, 9);
    mesh.addElement({f1_front, f2_front, f2_bottom, f2_top, f2_left, f2_right});

    // 2. 组装
    VEMSolver solver(mesh);
    Material mat; mat.E = 1000; mat.nu = 0.25;
    solver.assemble(mat);

    // 3. 设置边界条件 (Linear Patch Test)
    // 策略：将所有边界节点固定为精确解，
    // 检查内部节点 (如果有的话) 是否符合精确解。
    // 在这个 2-element 网格中，所有节点都在边界上... 
    // 这对于 Patch Test 有点作弊。
    // 为了真正测试，我们需要至少一个内部节点。
    // 让我们把上面的节点 4,5,6,7 (中间面) 视为 "自由" 节点，
    // 只固定左侧 (x=0) 和 右侧 (x=2) 的面。
    // 这样中间面的位移必须由 VEM 计算出来。
    
    std::map<int, double> dirichlet_bc;
    std::vector<int> fixed_nodes = {0,1,2,3, 8,9,10,11}; // Left and Right faces
    
    // 对于固定节点，施加精确解
    for (int nid : fixed_nodes) {
        Eigen::Vector3d u_exact = exactSolution(mesh.getNodes().col(nid));
        dirichlet_bc[3*nid + 0] = u_exact(0);
        dirichlet_bc[3*nid + 1] = u_exact(1);
        dirichlet_bc[3*nid + 2] = u_exact(2);
    }
    
    Eigen::VectorXd f_ext = Eigen::VectorXd::Zero(solver.getNumDofs());
    // Patch Test 无体力

    // 4. 求解
    Eigen::VectorXd u_sol = solver.solve(dirichlet_bc, f_ext);
    
    // 5. 验证误差 (Error Analysis)
    double max_error = 0.0;
    
    // 检查所有节点 (包括自由节点 4,5,6,7)
    for (int i = 0; i < mesh.getNodes().cols(); ++i) {
        Eigen::Vector3d u_num = u_sol.segment<3>(3*i);
        Eigen::Vector3d u_ex = exactSolution(mesh.getNodes().col(i));
        
        double err = (u_num - u_ex).norm();
        if (err > max_error) max_error = err;
        
        // Debug output for internal nodes
        if (i >= 4 && i <= 7) {
            // std::cout << "Node " << i << " error: " << err << std::endl;
        }
    }
    
    std::cout << "Max Displacement Error: " << max_error << std::endl;
    
    if (max_error < 1e-10) {
        std::cout << "[PASSED] Linear Patch Test passed!" << std::endl;
        return 0;
    } else {
        std::cerr << "[FAILED] Linear Patch Test failed!" << std::endl;
        return 1;
    }
}