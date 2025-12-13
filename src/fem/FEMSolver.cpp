#include "FEMSolver.hpp"
#include "LinearTet.hpp"
#include <iostream>
#include <set>
#include <vector>

namespace vem::fem {

FEMSolver::FEMSolver(const VEMMesh& mesh) : mesh_(mesh) {
    n_dofs_ = mesh_.getNumNodes() * 3;
}

void FEMSolver::assemble(const Material& mat) {
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(n_dofs_ * 20); // FEM 稀疏度比 VEM 低

    int n_elems = mesh_.getNumElements();
    
    for (int i = 0; i < n_elems; ++i) {
        // 1. 获取单元节点
        // VEMMesh 存储的是面。我们需要收集所有唯一节点。
        const auto& elem = mesh_.getElement(i);
        std::set<int> node_set;
        for (int fid : elem.face_indices) {
            const auto* face = mesh_.getFace(fid);
            for (int nid : face->node_indices) node_set.insert(nid);
        }
        
        // 检查是否为四面体
        if (node_set.size() != 4) {
            std::cerr << "[FEMSolver] Warning: Element " << i << " has " << node_set.size() 
                      << " nodes. Linear Tet FEM expects exactly 4. Skipping." << std::endl;
            continue;
        }

        std::vector<int> nodes(node_set.begin(), node_set.end());
        
        // 2. 准备坐标矩阵
        Eigen::MatrixXd elem_coords(3, 4);
        for (int k = 0; k < 4; ++k) {
            elem_coords.col(k) = mesh_.getNodes().col(nodes[k]);
        }

        // 3. 计算局部刚度矩阵
        Eigen::MatrixXd K_local = LinearTet::computeStiffness(elem_coords, mat);

        // 4. 组装到全局
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                // 3x3 block
                for (int di = 0; di < 3; ++di) {
                    for (int dj = 0; dj < 3; ++dj) {
                        int global_row = 3 * nodes[r] + di;
                        int global_col = 3 * nodes[c] + dj;
                        double val = K_local(3*r + di, 3*c + dj);
                        
                        if (std::abs(val) > 1e-16) {
                            triplets.emplace_back(global_row, global_col, val);
                        }
                    }
                }
            }
        }
    }

    K_global_.resize(n_dofs_, n_dofs_);
    K_global_.setFromTriplets(triplets.begin(), triplets.end());
}

Eigen::VectorXd FEMSolver::solve(const std::map<int, double>& dirichlet_bc, 
                                 const Eigen::VectorXd& f_ext) {
    // 这里为了偷懒且保证一致性，我建议完全复制 VEMSolver::solve 的逻辑
    // 或者将 Solver 逻辑抽象出来。
    // 现在直接复制 VEMSolver::solve 的代码即可，因为它只依赖 K_global_ 和 dof 索引
    // ... (请复制 VEMSolver::solve 实现，这里省略以节省篇幅，逻辑完全一样) ...
    // 如果你愿意，我可以把代码写全。
    
    if (K_global_.rows() == 0) return Eigen::VectorXd();

    // 1. Identify DOFs
    std::vector<int> free_dofs;
    std::vector<int> global_to_free(n_dofs_, -1);
    int n_fixed = 0;
    for(int i=0; i<n_dofs_; ++i) {
        if(dirichlet_bc.count(i)) n_fixed++;
        else {
            global_to_free[i] = free_dofs.size();
            free_dofs.push_back(i);
        }
    }
    int n_free = free_dofs.size();
    
    // 2. RHS
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_free);
    for(int i=0; i<n_free; ++i) if(free_dofs[i] < f_ext.size()) rhs(i) = f_ext(free_dofs[i]);

    // 3. Build K_FF and update RHS
    std::vector<Eigen::Triplet<double>> triplets_FF;
    for (int k=0; k<K_global_.outerSize(); ++k) {
        int c = k;
        bool col_fixed = (global_to_free[c] == -1);
        double u_val = col_fixed ? dirichlet_bc.at(c) : 0.0;
        
        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global_, k); it; ++it) {
            int r = it.row();
            if (global_to_free[r] != -1) { // Row is free
                if (col_fixed) rhs(global_to_free[r]) -= it.value() * u_val;
                else triplets_FF.emplace_back(global_to_free[r], global_to_free[c], it.value());
            }
        }
    }
    
    Eigen::SparseMatrix<double> K_FF(n_free, n_free);
    K_FF.setFromTriplets(triplets_FF.begin(), triplets_FF.end());
    
    // 4. Solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_FF);
    if(solver.info()!=Eigen::Success) return Eigen::VectorXd();
    Eigen::VectorXd u_F = solver.solve(rhs);
    
    Eigen::VectorXd u(n_dofs_);
    for(int i=0; i<n_free; ++i) u(free_dofs[i]) = u_F(i);
    for(auto const& [dof, val] : dirichlet_bc) u(dof) = val;
    return u;
}

} // namespace vem::fem