#pragma once
#include <Eigen/Sparse>
#include <map>
#include "../VEMMesh.hpp"
#include "../VEMElement.hpp" // For Material struct

namespace vem::fem {

class FEMSolver {
private:
    const VEMMesh& mesh_;
    Eigen::SparseMatrix<double> K_global_;
    int n_dofs_;

public:
    FEMSolver(const VEMMesh& mesh);

    // 组装 FEM 矩阵 (假设所有单元都是四面体)
    void assemble(const Material& mat);

    // 求解 (复用逻辑，或者 copy-paste 简化)
    Eigen::VectorXd solve(const std::map<int, double>& dirichlet_bc, 
                          const Eigen::VectorXd& f_ext);

    const Eigen::SparseMatrix<double>& getK() const { return K_global_; }
};

} // namespace vem::fem