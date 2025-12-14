#include "VEMSolver.hpp"
#include <iostream>
#include <set>
#include <vector>

VEMSolver::VEMSolver(const VEMMesh& mesh)
    : mesh_(mesh)
{
    n_dofs_ = mesh_.getNodes().cols() * 3; // 3 DOFs per node
}

void VEMSolver::assemble(const Material& mat)
{
    // 1. 准备 Triplets 用于构建稀疏矩阵
    // Triplet 格式: (row, col, value)
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;

    // 预估非零元素个数以加速内存分配 (经验值：每个节点关联约 20-50 个邻居?)
    // VEM 单元较密，保守估计
    triplets.reserve(n_dofs_ * 50);

    int n_elements = 0;
    // 假设已添加 getNumElements
    int num_elems = mesh_.getNumElements();

    for (int i = 0; i < num_elems; ++i) {
        if (!mesh_.isElementActive(i)) {
            continue;
        }
        const auto& elem = mesh_.getElement(i);

        // 2. 计算单元刚度矩阵
        VEMElement vem_elem(mesh_, elem);
        Eigen::MatrixXd K_el = vem_elem.computeStiffness(mat);

        std::set<int> unique_node_ids;
        for (int fid : elem.face_indices) {
            const auto* face = mesh_.getFace(fid);
            for (int nid : face->node_indices) {
                unique_node_ids.insert(nid);
            }
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
    const Eigen::VectorXd& f_ext)
{
    if (K_global_.rows() == 0) {
        std::cerr << "Error: Stiffness matrix not assembled!" << std::endl;
        return Eigen::VectorXd();
    }

    // 1. Identify Free and Fixed DOFs
    std::vector<int> free_dofs;
    std::vector<int> fixed_dofs;
    free_dofs.reserve(n_dofs_);

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
    if (n_free == 0) {
        // All DOFs fixed
        Eigen::VectorXd u(n_dofs_);
        for (auto const& [dof, val] : dirichlet_bc)
            u(dof) = val;
        return u;
    }

    // 2. Build Reduced System: K_FF * u_F = RHS
    // RHS = f_free - K_FD * u_D

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_free);

    // Add external forces
    for (int i = 0; i < n_free; ++i) {
        if (free_dofs[i] < f_ext.size())
            rhs(i) = f_ext(free_dofs[i]);
    }

    // Prepare Triplets for K_FF
    std::vector<Eigen::Triplet<double>> triplets_FF;
    triplets_FF.reserve(K_global_.nonZeros()); // Conservative estimate

    // Iterate over K_global to fill K_FF and update RHS
    // K_global is column-major by default
    for (int k = 0; k < K_global_.outerSize(); ++k) {
        int c = k; // Global Column Index

        bool col_is_fixed = (global_to_free[c] == -1);
        double u_col_val = 0.0;
        if (col_is_fixed) {
            u_col_val = dirichlet_bc.at(c);
        }

        for (Eigen::SparseMatrix<double>::InnerIterator it(K_global_, k); it; ++it) {
            int r = it.row(); // Global Row Index
            double val = it.value();

            bool row_is_free = (global_to_free[r] != -1);

            if (row_is_free) {
                int free_r = global_to_free[r];

                if (col_is_fixed) {
                    // Term belongs to K_FD * u_D
                    // Move to RHS: rhs[free_r] -= K_rc * u_c
                    rhs(free_r) -= val * u_col_val;
                } else {
                    // Term belongs to K_FF
                    int free_c = global_to_free[c];
                    triplets_FF.push_back(Eigen::Triplet<double>(free_r, free_c, val));
                }
            }
            // If row is fixed, we don't care (it's the reaction force equation)
        }
    }

    Eigen::SparseMatrix<double> K_FF(n_free, n_free);
    K_FF.setFromTriplets(triplets_FF.begin(), triplets_FF.end());

    std::cout << "add boundary condition done.\n";
    // 3. Solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(K_FF);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver factorization failed! Matrix might be singular or indefinite." << std::endl;
        return Eigen::VectorXd();
    }

    Eigen::VectorXd u_F = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Solver solution failed!" << std::endl;
        return Eigen::VectorXd();
    }

    // 4. Assemble full solution
    Eigen::VectorXd u = Eigen::VectorXd::Zero(n_dofs_);
    for (int i = 0; i < n_free; ++i)
        u(free_dofs[i]) = u_F(i);
    for (size_t i = 0; i < fixed_dofs.size(); ++i)
        u(fixed_dofs[i]) = dirichlet_bc.at(fixed_dofs[i]);

    return u;
}
