#include <iostream>
#include <vector>
#include <map>
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"

// 线性解析解 u(x,y,z)
// u = 0.1 + 0.1x
// v = 0.0
// w = 0.0
// 对应的应变 eps_x = 0.1, 其他 0. 
// 应力 sigma_x = (lambda+2mu)*0.1 + ...
// 这是一个非常简单的拉伸测试。
Eigen::Vector3d exactSolution(const Eigen::Vector3d& p) {
    return Eigen::Vector3d(0.1 * p.x(), 0.0, 0.0);
}

int main() {
    std::cout << "=== Single Element Linear Patch Test ===" << std::endl;

    // 1. 单个单位立方体
    VEMMesh mesh;
    Eigen::MatrixXd V(3, 8);
    V << 0, 1, 1, 0, 0, 1, 1, 0,
         0, 0, 1, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 1, 1, 1;
    mesh.setNodes(V);

    auto addFace = [&](int n1, int n2, int n3, int n4) {
        return mesh.addPolygonFace({n1, n2, n3, n4});
    };

    // 标准 6 面
    std::vector<int> faces;
    faces.push_back(addFace(0, 3, 2, 1)); // z=0
    faces.push_back(addFace(4, 5, 6, 7)); // z=1
    faces.push_back(addFace(0, 1, 5, 4)); // y=0
    faces.push_back(addFace(1, 2, 6, 5)); // x=1
    faces.push_back(addFace(2, 3, 7, 6)); // y=1
    faces.push_back(addFace(3, 0, 4, 7)); // x=0
    
    mesh.addElement(faces);

    // 2. 求解
    VEMSolver solver(mesh);
    Material mat; mat.E = 1000; mat.nu = 0.25;
    solver.assemble(mat);

    // 3. 边界条件
    // 固定所有 8 个节点为解析解 (Patch Test 的定义：边界位移固定，看内部是否一致)
    // 但单单元没有内部节点... 这变成了验证 K * u_exact = f_nodal?
    // 不，我们可以少固定几个点，让其他点自由，看能否算出正确位置。
    // 策略：固定 x=0 面 (nodes 0,3,4,7) 和 x=1 面 (nodes 1,2,5,6)。
    // 实际上这就是固定了所有节点...
    
    // 为了真正测试求解器，我们只固定能消除刚体位移的最少节点，
    // 或者固定 x=0 面，并在 x=1 面施加力？不，Patch Test 是 Dirichlet 问题。
    
    // 让我们做 "Displacement Recovery" 测试：
    // 如果我们固定 4 个节点（足够消除刚体位移），其他 4 个节点应该能算出来。
    // 固定 x=0 面 (Nodes 0,3,4,7) -> u=0
    // 固定 x=1 面 (Nodes 1,2,5,6) -> u=0.1
    // 等等，这就是简单的拉伸。
    
    // 让我们只固定 x=0 面完全固定。x=1 面施加 Neumann? 
    // 不，Patch Test 必须是全 Dirichlet 或混合。
    
    // 让我们尝试：固定 7 个节点，解第 8 个节点。
    // 这是一个极简的 Patch Test。
    std::map<int, double> dirichlet_bc;
    for (int i = 0; i < 8; ++i) {
        if (i == 6) continue; // 让节点 6 (1,1,1) 自由
        Eigen::Vector3d u_ex = exactSolution(mesh.getNodes().col(i));
        dirichlet_bc[3*i+0] = u_ex(0);
        dirichlet_bc[3*i+1] = u_ex(1);
        dirichlet_bc[3*i+2] = u_ex(2);
    }
    
    // 没有外力
    Eigen::VectorXd f_ext = Eigen::VectorXd::Zero(24);
    
    Eigen::VectorXd u_sol = solver.solve(dirichlet_bc, f_ext);
    
    // 检查节点 6 的结果
    Eigen::Vector3d u6_num = u_sol.segment<3>(3*6);
    Eigen::Vector3d u6_ex = exactSolution(mesh.getNodes().col(6));
    
    std::cout << "Node 6 Exact: " << u6_ex.transpose() << std::endl;
    std::cout << "Node 6 Num  : " << u6_num.transpose() << std::endl;
    
    double error = (u6_num - u6_ex).norm();
    std::cout << "Error: " << error << std::endl;
    
    if (error < 1e-10) {
        std::cout << "[PASSED] Single Element Patch Test" << std::endl;
        return 0;
    } else {
        std::cout << "[FAILED] Single Element Patch Test" << std::endl;
        return 1;
    }
}