#include "VEMElement.hpp"
#include <iostream>
#include <set>

Eigen::Matrix6d Material::getC() const {
    double l = lambda();
    double m = mu();
    
    Eigen::Matrix6d C = Eigen::Matrix6d::Zero();
    // 对角线
    C(0,0) = C(1,1) = C(2,2) = l + 2*m;
    C(3,3) = C(4,4) = C(5,5) = m;
    
    // 非对角线 (正应力耦合)
    C(0,1) = C(1,0) = C(0,2) = C(2,0) = C(1,2) = C(2,1) = l;
    
    return C;
}

VEMElement::VEMElement(const VEMMesh& mesh, const PolyhedronElement& elem)
    : mesh_(mesh), elem_(elem) {
    // 初始化缩放参数
    // scaling_center_ 取单元质心
    scaling_center_ = elem.centroid;
    
    // scaling_h_ 取单元体积的立方根，或者直径
    // 简单起见，使用 volume^(1/3)
    scaling_h_ = std::pow(elem.volume, 1.0/3.0);
}

Eigen::MatrixXd VEMElement::computeG(const Material& mat) const {
    int n_monos = getNumMonomials();
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(n_monos, n_monos);
    Eigen::Matrix6d C = mat.getC();
    double volume = elem_.volume;

    // 定义 helper lambda：计算第 i 个单项式基函数的常数应变 (Voigt notation 6x1)
    // Basis order (total 12):
    // 0-2:  [1,0,0], [0,1,0], [0,0,1]
    // 3-5:  [x~,0,0], [0,x~,0], [0,0,x~]  (x~ = (x-xc)/h)
    // 6-8:  [y~,0,0], [0,y~,0], [0,0,y~]
    // 9-11: [z~,0,0], [0,z~,0], [0,0,z~]
    auto getStrain = [&](int i) -> Eigen::Vector6d {
        Eigen::Vector6d eps = Eigen::Vector6d::Zero();
        double inv_h = 1.0 / scaling_h_;

        // 常数项 (0-2) 的梯度为0，应变也为0
        if (i < 3) return eps;

        // 线性项
        int linear_idx = i - 3; // 0..8
        int dir = linear_idx % 3; // 0=x, 1=y, 2=z (位移分量方向)
        int var = linear_idx / 3; // 0=x, 1=y, 2=z (变量依赖)

        // u = [.., .., ..]
        // du_i / dx_j 
        // strain_xx = du/dx, strain_yy = dv/dy, ...
        // strain_xy = 0.5 * (du/dy + dv/dx)
        
        // 我们的基函数是单分量非零。例如 [x~, 0, 0] -> u=(x-xc)/h, v=0, w=0
        // 只有 du/dx = 1/h 非零。
        
        // 填充 3x3 梯度矩阵 grad(u)
        Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();
        grad(dir, var) = inv_h; // du_dir / d_var

        // 转换为 Voigt 应变
        // eps = [xx, yy, zz, 2yz, 2xz, 2xy]
        eps(0) = grad(0,0); // xx
        eps(1) = grad(1,1); // yy
        eps(2) = grad(2,2); // zz
        eps(3) = grad(1,2) + grad(2,1); // 2yz
        eps(4) = grad(0,2) + grad(2,0); // 2xz
        eps(5) = grad(0,1) + grad(1,0); // 2xy
        
        return eps;
    };

    // 计算 G_ab = Vol * eps_a^T * C * eps_b
    // 注意：对于 k=1，P1 基函数的应变是常数。所以积分直接提出来变成体积。
    for (int i = 0; i < n_monos; ++i) {
        Eigen::Vector6d eps_i = getStrain(i);
        // 如果应变为0（刚体平移），G 对应的行/列全为 0
        if (eps_i.isZero()) continue;

        for (int j = i; j < n_monos; ++j) {
            Eigen::Vector6d eps_j = getStrain(j);
            if (eps_j.isZero()) continue;

            double val = volume * eps_i.dot(C * eps_j);
            G(i, j) = val;
            G(j, i) = val;
        }
    }

    return G;
}

Eigen::MatrixXd VEMElement::computeB(const Material& mat) const {
    int n_monos = getNumMonomials(); // 12
    
    // 1. 收集单元的所有唯一顶点，并建立映射
    std::set<int> unique_node_ids;
    for (int fid : elem_.face_indices) {
        const auto* face = mesh_.getFace(fid);
        for (int nid : face->node_indices) {
            unique_node_ids.insert(nid);
        }
    }
    
    // 建立 Global ID -> Local ID (0..N-1) 的映射
    // 同时也保存一个 vector 用于反查
    std::vector<int> local_nodes(unique_node_ids.begin(), unique_node_ids.end());
    std::map<int, int> global_to_local;
    for (size_t i = 0; i < local_nodes.size(); ++i) {
        global_to_local[local_nodes[i]] = i;
    }
    
    int n_nodes = local_nodes.size();
    int n_dofs = 3 * n_nodes;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(n_monos, n_dofs);
    
    Eigen::Matrix6d C = mat.getC();

    // 2. 遍历所有面进行积分
    for (int fid : elem_.face_indices) {
        const auto* face = mesh_.getFace(fid);
        
        // 计算该面上形状函数的积分权重: int_f phi_k dS
        auto shape_integrals = face->computeShapeFuncIntegrals(mesh_.getNodes());
        
        // 面的法向 (假设 computeProps 返回的是外法向，或者我们需要确保它是外法向)
        // 注意：VEMFace 中存储的法向是局部定义的。对于凸多面体，通常需要检查法向是否指向单元外部。
        // 一个简单的方法是检查 (FaceCentroid - ElemCentroid) . Normal > 0
        GeometricProps f_props = face->computeProps(mesh_.getNodes());
        Eigen::Vector3d outward_normal = f_props.normal;
        if ((f_props.centroid - elem_.centroid).dot(outward_normal) < 0) {
            outward_normal = -outward_normal;
        }

        // 3. 遍历所有单项式 m_alpha (Rows of B)
        // 我们复用 computeG 中的 lambda (需要把它提出来或者复制一遍)
        // 为了简单，建议把 getStrain 变成 VEMElement 的私有辅助函数
        // 这里我先复制一遍逻辑
        auto getStrain = [&](int i) -> Eigen::Vector6d {
             // ... (复制之前的 getStrain 代码) ...
             // 请务必保证 scaling_h_ 等参数一致
             Eigen::Vector6d eps = Eigen::Vector6d::Zero();
             double inv_h = 1.0 / scaling_h_;
             if (i < 3) return eps;
             int linear_idx = i - 3;
             int dir = linear_idx % 3;
             int var = linear_idx / 3;
             Eigen::Matrix3d grad = Eigen::Matrix3d::Zero();
             grad(dir, var) = inv_h; 
             eps(0) = grad(0,0); eps(1) = grad(1,1); eps(2) = grad(2,2);
             eps(3) = grad(1,2) + grad(2,1); eps(4) = grad(0,2) + grad(2,0); eps(5) = grad(0,1) + grad(1,0);
             return eps;
        };

        for (int alpha = 0; alpha < n_monos; ++alpha) {
            Eigen::Vector6d eps_alpha = getStrain(alpha);
            
            // 如果应变为0 (刚体平移)，则 traction = 0，B 对应行也是 0
            if (eps_alpha.isZero()) continue;

            // 计算应力 sigma = C * eps
            Eigen::Vector6d sigma = C * eps_alpha;
            
            // 计算面上的牵引力 t = sigma * n (Voigt notation product)
            // t_x = sig_xx * nx + sig_xy * ny + sig_xz * nz
            // t_y = sig_yx * nx + sig_yy * ny + sig_yz * nz
            // t_z = sig_zx * nx + sig_zy * ny + sig_zz * nz
            
            Eigen::Vector3d traction;
            double nx = outward_normal.x();
            double ny = outward_normal.y();
            double nz = outward_normal.z();
            
            traction.x() = sigma(0)*nx + sigma(5)*ny + sigma(4)*nz;
            traction.y() = sigma(5)*nx + sigma(1)*ny + sigma(3)*nz;
            traction.z() = sigma(4)*nx + sigma(3)*ny + sigma(2)*nz;

            // 4. 累加到矩阵 B
            // B(alpha, dof) += int_f (t . v) dS
            // 由于 t 是常数，int_f t . v dS = t . (sum_k v_k * weight_k)
            // 所以对于每个节点 k，贡献是 t_dir * weight_k
            
            for (auto const& [node_id, weight] : shape_integrals) {
                int local_id = global_to_local[node_id];
                
                // DoF 排列: [node0_x, node0_y, node0_z, node1_x, ...]
                B(alpha, 3*local_id + 0) += traction.x() * weight;
                B(alpha, 3*local_id + 1) += traction.y() * weight;
                B(alpha, 3*local_id + 2) += traction.z() * weight;
            }
        }
    }

    return B;
}

Eigen::MatrixXd VEMElement::computeD() const {
    int n_monos = getNumMonomials();
    
    // 重新收集节点 (复用 computeB 中的逻辑，或者将其提取为辅助函数)
    // 这里为了独立性，再写一遍
    std::set<int> unique_node_ids;
    for (int fid : elem_.face_indices) {
        const auto* face = mesh_.getFace(fid);
        for (int nid : face->node_indices) unique_node_ids.insert(nid);
    }
    std::vector<int> local_nodes(unique_node_ids.begin(), unique_node_ids.end());
    
    int n_dofs = 3 * local_nodes.size();
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n_dofs, n_monos);
    
    double inv_h = 1.0 / scaling_h_;
    Eigen::Vector3d xc = scaling_center_;
    
    // 遍历所有节点
    for (size_t i = 0; i < local_nodes.size(); ++i) {
        Eigen::Vector3d pos = mesh_.getNodes().col(local_nodes[i]);
        Eigen::Vector3d diff = (pos - xc) * inv_h; // Normalized position
        
        // 遍历所有单项式
        // 0-2: Constant [1,0,0]...
        // 3-11: Linear
        for (int alpha = 0; alpha < n_monos; ++alpha) {
            Eigen::Vector3d val = Eigen::Vector3d::Zero();
            
            if (alpha < 3) {
                // Constant bases: 0->x, 1->y, 2->z
                val(alpha) = 1.0;
            } else {
                // Linear bases
                int linear_idx = alpha - 3;
                int dir = linear_idx % 3; // Component direction
                int var = linear_idx / 3; // Variable dependency
                
                // e.g. alpha corresponding to u = x_tilde * e_x
                // dir = 0 (x), var = 0 (x)
                val(dir) = diff(var);
            }
            
            // Fill D matrix
            // DOF layout: [node0_x, node0_y, node0_z, node1_x, ...]
            D(3*i + 0, alpha) = val(0);
            D(3*i + 1, alpha) = val(1);
            D(3*i + 2, alpha) = val(2);
        }
    }
    
    return D;
}

Eigen::MatrixXd VEMElement::computeStiffness(const Material& mat) const {
    // 1. Compute Matrices
    Eigen::MatrixXd G = computeG(mat);
    Eigen::MatrixXd B = computeB(mat);
    Eigen::MatrixXd D = computeD();
    
    // 2. Consistency Term: K_c = B^T * G^+ * B
    // G is singular (rank deficiency 6), so we use Pseudo-Inverse
    Eigen::MatrixXd G_pinv = G.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::MatrixXd K_c = B.transpose() * G_pinv * B;
    
    // 3. Stabilization Term: K_s
    // We use a least-squares projection stabilization
    // P_ls = (D^T D)^-1 D^T
    // Stabilizer = (I - D * P_ls)^T * (I - D * P_ls)
    // This penalizes the difference between u and its best-fit polynomial P(u)
    
    // D^T * D is small (12x12), easy to invert
    Eigen::MatrixXd DtD = D.transpose() * D;
    Eigen::MatrixXd P_ls = DtD.ldlt().solve(D.transpose()); // Robust solve
    
    int n_dofs = D.rows();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n_dofs, n_dofs);
    Eigen::MatrixXd Proj_defect = I - D * P_ls;
    
    Eigen::MatrixXd K_s_raw = Proj_defect.transpose() * Proj_defect;
    
    // Calculate stabilization parameter alpha
    // Typically trace(K_c) / N_dofs or similar scaling
    // Use standard trace scaling
    double trace_Kc = K_c.trace();
    // Avoid division by zero if K_c is empty or zero (though unlikely)
    double alpha = (trace_Kc > 1e-12) ? (trace_Kc / n_dofs) : 1.0; 
    
    // Usually stability scaling factor is around 0.5 to 1.0 times the mean diagonal stiffness
    // Let's use 1.0 * mean_diag
    // Gain 2014 suggests stability parameter based on material trace, but trace method is robust.
    
    Eigen::MatrixXd K_s = alpha * K_s_raw;
    
    return K_c + K_s;
}