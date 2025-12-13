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