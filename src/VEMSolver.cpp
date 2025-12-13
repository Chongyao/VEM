#include "VEMSolver.hpp"
#include <iostream>
#include <set>
#include <vector>

VEMSolver::VEMSolver(const VEMMesh& mesh) : mesh_(mesh) {
    n_dofs_ = mesh_.getNodes().cols() * 3; // 3 DOFs per node
}

void VEMSolver::assemble(const Material& mat) {
    // 1. 准备 Triplets 用于构建稀疏矩阵
    // Triplet 格式: (row, col, value)
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    
    // 预估非零元素个数以加速内存分配 (经验值：每个节点关联约 20-50 个邻居?)
    // VEM 单元较密，保守估计
    triplets.reserve(n_dofs_ * 50);

    int n_elements = 0; // 可以添加接口获取单元数量，目前只能通过 exportToVTK 看到，或者我们在 Mesh 中加个 getNumElements
    // 修正：VEMMesh 目前没有直接暴露 elements_ 向量大小的接口。
    // 我们需要在 VEMMesh.hpp 中添加 int getNumElements() const { return elements_.size(); }
    // 假设您会去添加这个接口，或者我们现在先用 getElement(i) 尝试直到失败？
    // 为了代码稳健，请先去 VEMMesh.hpp 添加 getNumElements()。
    
    // 假设已添加 getNumElements
    int num_elems = mesh_.getNumElements();

    for (int i = 0; i < num_elems; ++i) {
        const auto& elem = mesh_.getElement(i);
        
        // 2. 计算单元刚度矩阵
        VEMElement vem_elem(mesh_, elem);
        Eigen::MatrixXd K_el = vem_elem.computeStiffness(mat);
        
        // 3. 获取单元节点的全局索引
        // 我们需要复用 "收集节点" 的逻辑。
        // 为了性能，最好在 VEMElement 中暴露 "getLocalToGlobalNodeMap" 
        // 但现在我们重新收集一遍 (逻辑必须与 computeStiffness 内部完全一致!)
        // 
        // 注意：VEMElement::computeD/B 内部是按 unique_node_ids (std::set 排序) 来排列 DOF 的。
        // 我们必须保证这里也是一样的顺序。
        
        std::set<int> unique_node_ids;
        for (int fid : elem.face_indices) {
            const auto* face = mesh_.getFace(fid);
            for (int nid : face->node_indices) unique_node_ids.insert(nid);
        }
        std::vector<int> local_nodes(unique_node_ids.begin(), unique_node_ids.end());
        
        int n_local_nodes = local_nodes.size();
        
        // 4. 填充 Triplets
        for (int r = 0; r < n_local_nodes; ++r) {
            int global_node_r = local_nodes[r];
            
            for (int c = 0; c < n_local_nodes; ++c) {
                int global_node_c = local_nodes[c];
                
                // 3x3 block for each node pair
                for (int di = 0; di < 3; ++di) { // direction i (row)
                    for (int dj = 0; dj < 3; ++dj) { // direction j (col)
                        
                        int local_row = 3 * r + di;
                        int local_col = 3 * c + dj;
                        
                        int global_row = 3 * global_node_r + di;
                        int global_col = 3 * global_node_c + dj;
                        
                        double val = K_el(local_row, local_col);
                        
                        // 仅添加非零值 (可选，但推荐)
                        if (std::abs(val) > 1e-16) {
                            triplets.push_back(T(global_row, global_col, val));
                        }
                    }
                }
            }
        }
    }

    // 5. 从 Triplets 构建稀疏矩阵
    K_global_.resize(n_dofs_, n_dofs_);
    K_global_.setFromTriplets(triplets.begin(), triplets.end());
}


Eigen::VectorXd VEMSolver::solve(const std::map<int, double>& dirichlet_bc, 
                                 const Eigen::VectorXd& f_ext) {
    if (K_global_.rows() == 0) {
        std::cerr << "Error: Stiffness matrix not assembled!" << std::endl;
        return Eigen::VectorXd();
    }

    // 1. 识别自由 (Free) 和固定 (Dirichlet) 的自由度索引
    std::vector<int> free_dofs;
    std::vector<int> fixed_dofs;
    free_dofs.reserve(n_dofs_);
    
    // 我们需要一个快速查找表来构建映射
    std::vector<int> global_to_free(n_dofs_, -1);

    for (int i = 0; i < n_dofs_; ++i) {
        if (dirichlet_bc.find(i) != dirichlet_bc.end()) {
            fixed_dofs.push_back(i);
        } else {
            global_to_free[i] = free_dofs.size();
            free_dofs.push_back(i);
        }
    }

    int n_free = free_dofs.size();
    
    // 2. 构建缩减系统 K_FF * u_F = rhs_F
    // RHS = f_ext - K_FD * u_D
    
    // 构造 u_D 向量
    Eigen::VectorXd u_D(fixed_dofs.size());
    for (size_t i = 0; i < fixed_dofs.size(); ++i) {
        u_D(i) = dirichlet_bc.at(fixed_dofs[i]);
    }

    // 提取 K_FF 和 K_FD
    // 这一步对于大规模矩阵可能比较慢，可以使用 Eigen 的切片操作或重新组装
    // 简单起见，我们遍历 Triplet 或者直接遍历 K_global
    
    Eigen::SparseMatrix<double> K_FF(n_free, n_free);
    Eigen::SparseMatrix<double> K_FD(n_free, fixed_dofs.size());
    
    // 辅助: 从 global 提取子矩阵
    // 为了效率，我们使用 setFromTriplets
    std::vector<Eigen::Triplet<double>> triplets_FF;
    std::vector<Eigen::Triplet<double>> triplets_FD;
    
    // 遍历 K_global 的非零元素
    for (int k = 0; k < K_global_.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global_, k); it; ++it) {
            int r = it.row();
            int c = it.col(); // col index (here k)
            double val = it.value();
            
            bool r_is_free = (global_to_free[r] != -1);
            bool c_is_free = (global_to_free[c] != -1);
            
            if (r_is_free && c_is_free) {
                triplets_FF.push_back(Eigen::Triplet<double>(global_to_free[r], global_to_free[c], val));
            } else if (r_is_free && !c_is_free) {
                // K_FD 部分: 列是固定的
                // 我们需要找到 c 在 fixed_dofs 中的索引吗？
                // 其实不需要构建完整的 K_FD 矩阵，我们可以直接更新 RHS
                // RHS_i -= val * u_D_j
                // 这里为了清晰，我们在后面直接处理 RHS
            }
        }
    }
    
    K_FF.setFromTriplets(triplets_FF.begin(), triplets_FF.end());
    
    // 3. 计算 RHS
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_free);
    
    // 添加外力项 f_free
    for (int i = 0; i < n_free; ++i) {
        rhs(i) = f_ext(free_dofs[i]);
    }
    
    // 减去 K_FD * u_D
    // 我们可以直接遍历 K_global 来做这个操作，比构建 K_FD 更快
    for (int k = 0; k < K_global_.outerSize(); ++k) {
        int c = k; // global col
        bool c_is_fixed = (global_to_free[c] == -1);
        
        if (c_is_fixed) {
            double u_val = dirichlet_bc.at(c);
            
            for (Eigen::SparseMatrix<double>::InnerIterator it(K_global_, k); it; ++it) {
                int r = it.row();
                double val = it.value();
                
                if (global_to_free[r] != -1) {
                    // Row is free, Col is fixed
                    // rhs[free_row] -= K_rc * u_c
                    rhs(global_to_free[r]) -= val * u_val;
                }
            }
        }
    }

    // 4. 求解
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_FF);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver factorization failed!" << std::endl;
        return Eigen::VectorXd();
    }
    
    Eigen::VectorXd u_F = solver.solve(rhs);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver solution failed!" << std::endl;
        return Eigen::VectorXd();
    }

    // 5. 组装完整解向量
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n_dofs_);
    for (int i = 0; i < n_free; ++i) u(free_dofs[i]) = u_F(i);
    for (size_t i = 0; i < fixed_dofs.size(); ++i) u(fixed_dofs[i]) = u_D(i);
    
    return u;
}