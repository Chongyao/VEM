#include <iostream>
#include <set>
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

    // 3. 计算矩阵
    Eigen::MatrixXd G = vem_elem.computeG(mat);
    Eigen::MatrixXd B = vem_elem.computeB(mat); // New!
    
    std::cout << "Computed G size: " << G.rows() << "x" << G.cols() << std::endl;
    std::cout << "Computed B size: " << B.rows() << "x" << B.cols() << std::endl;

    // 获取节点的局部映射 (为了构建测试位移向量 u)
    // 我们需要知道 computeB 内部的节点顺序。
    // 由于 computeB 内部是按 unique_node_ids (std::set 自动排序) 排列的，
    // 这里我们用同样的逻辑重建映射。
    std::set<int> unique_node_ids;
    for (int fid : elem.face_indices) {
        for (int nid : mesh.getFace(fid)->node_indices) unique_node_ids.insert(nid);
    }
    std::vector<int> local_nodes(unique_node_ids.begin(), unique_node_ids.end());


    // ==========================================
    // 验证 1: B * u_RBM = 0 (刚体运动检验)
    // ==========================================
    // 构造一个刚体平移 u = [1, 0, 0]
    Eigen::VectorXd u_trans = Eigen::VectorXd::Zero(B.cols());
    for (size_t i = 0; i < local_nodes.size(); ++i) {
        u_trans(3*i + 0) = 1.0; 
    }
    
    // B * u 应该近似为 0
    double err_trans = (B * u_trans).norm();
    std::cout << "Error B * u_trans: " << err_trans << std::endl;
    checkCondition(err_trans < 1e-9, "B matrix passes Rigid Body Translation test.");


    // ==========================================
    // 验证 2: G 矩阵的一致性检验 (Polynomial Consistency)
    // B * u_beta = G_col_beta
    // 如果 u 对应第 beta 个单项式，那么 B*u 应该等于 G 的第 beta 列
    // ==========================================
    
    // 选取一个非零能量的单项式，例如 index 3: m = (x-xc)/h * e_x (对应 u=(x-xc)/h, v=0, w=0)
    // 这是一个线性拉伸模式
    int beta = 3; 
    Eigen::VectorXd u_poly = Eigen::VectorXd::Zero(B.cols());
    
    // 单元质心和缩放 (需要与 VEMElement 内部一致，这里手动计算一下)
    Eigen::Vector3d xc = elem.centroid;
    double h = std::pow(elem.volume, 1.0/3.0);
    
    // 构造节点位移场
    for (size_t i = 0; i < local_nodes.size(); ++i) {
        int nid = local_nodes[i];
        Eigen::Vector3d pos = mesh.getNodes().col(nid);
        
        // 单项式 3 定义: vector field = [(x-xc)/h, 0, 0]
        u_poly(3*i + 0) = (pos.x() - xc.x()) / h;
        u_poly(3*i + 1) = 0.0;
        u_poly(3*i + 2) = 0.0;
    }
    
    Eigen::VectorXd Bu = B * u_poly;
    Eigen::VectorXd G_col = G.col(beta);
    
    double err_poly = (Bu - G_col).norm();
    
    // 打印对比 (可选)
    // std::cout << "B*u: " << Bu.transpose() << std::endl;
    // std::cout << "G_col: " << G_col.transpose() << std::endl;
    
    std::cout << "Error B * u_poly - G_col: " << err_poly << std::endl;
    // 注意：这里的误差可能略大一点点，取决于数值积分精度，但对于线性函数应该是机器精度
    checkCondition(err_poly < 1e-9, "B matrix satisfies Polynomial Consistency (B*u = G*s).");

    return 0;
}