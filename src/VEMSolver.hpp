#pragma once
#include <Eigen/Sparse>
#include "VEMMesh.hpp"
#include "VEMElement.hpp"

class VEMSolver {
private:
    const VEMMesh& mesh_;
    
    // 全局刚度矩阵 (稀疏)
    Eigen::SparseMatrix<double> K_global_;
    
    // 自由度总数 (3 * NumNodes)
    int n_dofs_;

public:
    VEMSolver(const VEMMesh& mesh);

    // 组装全局刚度矩阵
    void assemble(const Material& mat);

    // 获取组装好的稀疏矩阵 (用于调试或后续求解)
    const Eigen::SparseMatrix<double>& getK() const { return K_global_; }
    
    // 获取总自由度数
    int getNumDofs() const { return n_dofs_; }

    // 新增：求解线性系统 K u = f
    // dirichlet_bc: map<global_dof_index, value>
    // f_ext: 外部力向量 (大小为 n_dofs)
    // 返回: 全局位移向量 u
    Eigen::VectorXd solve(const std::map<int, double>& dirichlet_bc, 
                          const Eigen::VectorXd& f_ext);

};