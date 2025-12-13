#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"

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
    // ==========================================
    // 1. 构建网格：两个并排的单位立方体
    // Cube 1: [0,1] x [0,1] x [0,1]
    // Cube 2: [1,2] x [0,1] x [0,1]
    // 共享面 x=1
    // ==========================================
    VEMMesh mesh;

    // 节点列表 (12个节点)
    // 0-3: x=0 面 (Left)
    // 4-7: x=1 面 (Middle, Shared)
    // 8-11: x=2 面 (Right)
    Eigen::MatrixXd V(3, 12);
    V << 0,0,0,0, 1,1,1,1, 2,2,2,2, // x
         0,1,1,0, 0,1,1,0, 0,1,1,0, // y
         0,0,1,1, 0,0,1,1, 0,0,1,1; // z
    mesh.setNodes(V);

    // 定义面 (需要小心顶点的顺序以保证法向朝外，或者依赖代码内部的自动修正)
    // 这里我们手动定义两个单元的面
    
    // 辅助lambda：添加一个四边形面
    auto addFace = [&](int n1, int n2, int n3, int n4) {
        return mesh.addPolygonFace({n1, n2, n3, n4});
    };

    // --- Element 1 Faces ---
    std::vector<int> faces1;
    faces1.push_back(addFace(0, 1, 2, 3)); // z=1 ? No, check coords. 
    // 让我们更仔细地定义拓扑
    // y-z 平面投影: 0(00), 1(10), 2(11), 3(01). 
    // x=0: 0-3-2-? 让我们用简单的顺序，VEMMesh会自动计算几何
    // 为了简单，我们只保证面的连接关系是对的
    
    // 重新定义面的节点索引以确保闭合
    // Cube 1 faces
    int f1_back   = addFace(0, 3, 2, 1); // x=0 (Normal -x)
    int f1_front  = addFace(4, 5, 6, 7); // x=1 (Normal +x) -> shared
    int f1_bottom = addFace(0, 1, 5, 4); // z=0
    int f1_top    = addFace(3, 7, 6, 2); // z=1
    int f1_left   = addFace(0, 4, 7, 3); // y=0
    int f1_right  = addFace(1, 2, 6, 5); // y=1
    
    mesh.addElement({f1_back, f1_front, f1_bottom, f1_top, f1_left, f1_right});

    // --- Element 2 Faces ---
    // 注意：Cube 2 共享 f1_front (x=1) 作为它的“背面”
    // 但在 VEM 中，通常每个单元引用自己的面列表。
    // 如果两个单元共享同一个面ID，VEMMesh 需要能处理（我们的实现是基于面的，所以共享面ID是完全合法的，也是必须的）
    // 共享面: f1_front (4,5,6,7). 对于 Cube 1 它是前，对于 Cube 2 它是后。
    
    int f2_front  = addFace(8, 9, 10, 11); // x=2
    int f2_bottom = addFace(4, 5, 9, 8);   // z=0
    int f2_top    = addFace(7, 11, 10, 6); // z=1
    int f2_left   = addFace(4, 8, 11, 7);  // y=0
    int f2_right  = addFace(5, 6, 10, 9);  // y=1
    
    // 关键：重用 f1_front
    mesh.addElement({f1_front, f2_front, f2_bottom, f2_top, f2_left, f2_right});

    std::cout << "Mesh created: " << mesh.getNodes().cols() << " nodes, " 
              << mesh.getNumElements() << " elements." << std::endl;

    // ==========================================
    // 2. 组装全局刚度矩阵
    // ==========================================
    VEMSolver solver(mesh);
    
    Material mat;
    mat.E = 1000.0;
    mat.nu = 0.3;

    solver.assemble(mat);
    
    const Eigen::SparseMatrix<double>& K_global = solver.getK();

    std::cout << "Global Stiffness Matrix: " << K_global.rows() << "x" << K_global.cols() 
              << ", Non-zeros: " << K_global.nonZeros() << std::endl;

    // 验证大小
    checkCondition(K_global.rows() == 36 && K_global.cols() == 36, "Global K size is correct (36x36).");

    // ==========================================
    // 3. 验证刚体模式 (Global Nullity)
    // 整个连通结构应该只有 6 个零特征值
    // ==========================================
    
    // 将稀疏矩阵转换为稠密矩阵进行特征值分解 (仅适用于小规模测试)
    Eigen::MatrixXd K_dense = Eigen::MatrixXd(K_global);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K_dense);
    Eigen::VectorXd ev = es.eigenvalues();

    int zero_ev = 0;
    int neg_ev = 0;
    for (int i = 0; i < ev.size(); ++i) {
        if (std::abs(ev(i)) < 1e-8) zero_ev++;
        else if (ev(i) < -1e-8) neg_ev++;
    }

    std::cout << "Global Eigenvalues (first 10): " << ev.head(10).transpose() << std::endl;

    checkCondition(neg_ev == 0, "Global K is positive semi-definite.");
    
    if (zero_ev == 6) {
        std::cout << "[PASSED] Global K has exactly 6 Rigid Body Modes (Elements are connected properly)." << std::endl;
    } else {
        std::cerr << "[FAILED] Expected 6 Rigid Body Modes, found " << zero_ev << "." << std::endl;
        if (zero_ev == 12) {
            std::cerr << "Hint: It looks like the elements are not connected (12 RBMs)." << std::endl;
        }
        return 1;
    }

    return 0;
}