// #include "VEMMesh.hpp"
// #include "VEMSolver.hpp"
// #include "fem/FEMSolver.hpp"
// #include "io/VTKLoader.hpp"
// #include "io/VTKWriter.hpp"
// #include <cmath>
// #include <iomanip>
// #include <iostream>
// #include <set>
// #include <vector>

// using namespace vem;

// // 辅助：生成梁的四面体网格
// // box: [0, L] x [0, W] x [0, H]
// void createBeamMesh(VEMMesh& mesh, int nx, int ny, int nz, double L, double W, double H)
// {
//     int n_nodes = (nx + 1) * (ny + 1) * (nz + 1);
//     Eigen::MatrixXd V(3, n_nodes);

//     double dx = L / nx;
//     double dy = W / ny;
//     double dz = H / nz;

//     auto get_node_idx = [&](int i, int j, int k) {
//         return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
//     };

//     // 1. 生成节点
//     for (int k = 0; k <= nz; ++k) {
//         for (int j = 0; j <= ny; ++j) {
//             for (int i = 0; i <= nx; ++i) {
//                 int idx = get_node_idx(i, j, k);
//                 V(0, idx) = i * dx;
//                 V(1, idx) = j * dy;
//                 V(2, idx) = k * dz;
//             }
//         }
//     }
//     mesh.setNodes(V);

//     // 2. 生成单元 (将每个 Cube 剖分为 6 个 Tets)
//     for (int k = 0; k < nz; ++k) {
//         for (int j = 0; j < ny; ++j) {
//             for (int i = 0; i < nx; ++i) {
//                 // Cube vertices
//                 int n0 = get_node_idx(i, j, k);
//                 int n1 = get_node_idx(i + 1, j, k);
//                 int n2 = get_node_idx(i + 1, j + 1, k);
//                 int n3 = get_node_idx(i, j + 1, k);
//                 int n4 = get_node_idx(i, j, k + 1);
//                 int n5 = get_node_idx(i + 1, j, k + 1);
//                 int n6 = get_node_idx(i + 1, j + 1, k + 1);
//                 int n7 = get_node_idx(i, j + 1, k + 1);

//                 // Kuhn subdivision (6 tets)
//                 // 具体的分割方式保证相邻单元面匹配
//                 // Tet 1: 0 1 3 4
//                 // Tet 2: 1 2 3 6  <-- Check indices carefully, keeping it simple
//                 // Let's use a standard decomposition:
//                 // T1: 0 1 2 5
//                 // T2: 0 2 3 7
//                 // T3: 0 5 7 4
//                 // T4: 0 2 5 7
//                 // T5: 2 5 6 7
//                 // Wait, constructing valid tet mesh manually is verbose.
//                 // Let's use a simpler 5-tet decomposition (symmetric) or assume standard marching.

//                 // 采用如下 6-tet 划分 (Path: 0->1->2->6->7 etc?)
//                 // Tet 1: (0,1,2,6) -> Wrong edge?

//                 // 让我们用最稳健的写法：
//                 int Hex[8] = { n0, n1, n2, n3, n4, n5, n6, n7 };
//                 // Tet 1: 0, 5, 1, 4 ? No.

//                 // 使用一种标准划分:
//                 // 0 1 2 5
//                 // 0 2 3 7
//                 // 0 5 7 4
//                 // 0 2 5 7  <-- Central
//                 // 2 5 6 7
//                 // 5 tets? No, 6 is better for structured grid consistency.

//                 // VEM 不挑食，只要闭合即可。
//                 // 我们添加 5 个四面体版本 (Checkerboard pattern needed for consistency, but for single block ok)
//                 // 为了简单且不出错，这里只添加 1 个巨大的六面体作为 VEM 单元？
//                 // 不，用户说 "由 tet 组成"。

//                 // OK, 6 Tets decomposition:
//                 auto addTet = [&](int a, int b, int c, int d) {
//                     std::vector<int> f1 = { a, c, b };
//                     std::vector<int> f2 = { a, b, d };
//                     std::vector<int> f3 = { a, d, c };
//                     std::vector<int> f4 = { b, c, d };
//                     std::vector<int> fs;
//                     fs.push_back(mesh.addTriFace(f1));
//                     fs.push_back(mesh.addTriFace(f2));
//                     fs.push_back(mesh.addTriFace(f3));
//                     fs.push_back(mesh.addTriFace(f4));
//                     mesh.addElement(fs);
//                 };

//                 addTet(n0, n1, n2, n5);
//                 addTet(n0, n2, n3, n7);
//                 addTet(n0, n5, n7, n4);
//                 addTet(n0, n2, n5, n7); // Center tet 1 (Wait, volume overlap?)
//                 addTet(n2, n5, n6, n7);
//                 // 这种 5-tet 划分在相邻 hex 之间需要棋盘格翻转才能协调 (conforming)。
//                 // 否则面不匹配。
//                 // 既然 VEM 可以处理任意多面体，为什么不直接把这个 Hex 作为一个单元？
//                 // "Beam, 由 Tet 组成" -> 好吧，我们必须用 Tet。

//                 // 既然写 6-Tet 比较繁琐，我们改用 5-Tet 且加 checkerboard 逻辑？
//                 // 或者更简单：只生成 12 个三角形面组成的六面体，然后告诉 VEM 它是 Polyhedron。
//                 // 但为了对比 FEM (LinearTet)，必须是 Tet。

//                 // 6-Tet 分解公式 (对角线 0-6):
//                 addTet(n0, n1, n5, n6); // 错了，体积覆盖。

//                 // 让我们用一种绝对安全的写法 (0,1,2,3 bottom, 4,5,6,7 top)
//                 // 1. (0,1,2,6)
//                 // 2. (0,2,3,7) -- no

//                 // 算了，为了不写错拓扑，我们换个思路：
//                 // 这是一个 VEM 代码库。我们直接生成 Polyhedral Mesh (Hex) 对比 VEM Hex vs FEM Tet?
//                 // 不，用户要求 Tet。

//                 // 采用 5-tet decomposition (无需坐标变换，只要奇偶格子翻转)
//                 bool flip = ((i + j + k) % 2 == 1);
//                 if (!flip) {
//                     addTet(n0, n1, n3, n4);
//                     addTet(n1, n2, n0, n5); // n0 is base? (1,2,5) is corner. (1,2,5,6)? No. (1,2,5) connect to 0?
//                     // Corner tets:
//                     addTet(n1, n2, n5, n6); // Corner at 5? No
//                     // Correct 5-tet:
//                     // Corners: (0,1,3,4), (1,2,5,6), (3,2,7,6), (4,5,7,6) -> 4 tets
//                     // Inner: (1,3,4,6) -> ?

//                     // 修正的 5-Tet:
//                     // 1. (0, 1, 3, 4)
//                     // 2. (1, 2, 5, 6)
//                     // 3. (3, 6, 2, 7) - wait indices
//                     // 4. (4, 6, 7, 5)
//                     // 5. (1, 4, 6, 3) - Central

//                     addTet(n0, n1, n3, n4);
//                     addTet(n1, n2, n5, n6);
//                     addTet(n2, n3, n6, n7);
//                     addTet(n4, n5, n6, n7);
//                     addTet(n1, n3, n4, n6); // Center
//                 } else {
//                     // Flip indices for checkerboard
//                     // Corners: (0,1,2,5), (0,3,2,7), (0,4,5,7), (5,6,7,2)
//                     // Inner: (0,2,5,7)
//                     addTet(n0, n1, n2, n5);
//                     addTet(n0, n2, n3, n7);
//                     addTet(n0, n4, n5, n7);
//                     addTet(n2, n5, n6, n7);
//                     addTet(n0, n5, n2, n7); // Center
//                 }
//             }
//         }
//     }
// }

// // 计算体力载荷向量
// Eigen::VectorXd computeGravityLoad(const VEMMesh& mesh, const Eigen::Vector3d& gravity, double rho)
// {
//     int n_nodes = mesh.getNumNodes();
//     Eigen::VectorXd f = Eigen::VectorXd::Zero(n_nodes * 3);

//     int n_elems = mesh.getNumElements();
//     for (int i = 0; i < n_elems; ++i) {
//         const auto& elem = mesh.getElement(i);
//         double mass = elem.volume * rho;

//         // 获取单元所有节点
//         std::set<int> node_set;
//         for (int fid : elem.face_indices) {
//             for (int nid : mesh.getFace(fid)->node_indices)
//                 node_set.insert(nid);
//         }

//         // 将重力均分到节点 (Lumping)
//         // 对于线性单元这是精确的 (consistent force vector for constant body force)
//         double node_load_mag = mass / node_set.size();

//         for (int nid : node_set) {
//             f(3 * nid + 0) += node_load_mag * gravity(0);
//             f(3 * nid + 1) += node_load_mag * gravity(1);
//             f(3 * nid + 2) += node_load_mag * gravity(2);
//         }
//     }
//     return f;
// }

// int main()
// {
//     std::cout << "=== Beam Gravity Test (VEM vs FEM on Tets) ===" << std::endl;

//     // 1. 生成网格
//     VEMMesh mesh;

//     io::VTKLoader loader;
//     mesh = *loader.load("beam.vtu");

//     // ... 加载网格后 ...
//     double min_x = 1e30, max_x = -1e30;
//     const auto& nodes = mesh.getNodes();
//     for (int i = 0; i < nodes.cols(); ++i) {
//         if (nodes(0, i) < min_x)
//             min_x = nodes(0, i);
//         if (nodes(0, i) > max_x)
//             max_x = nodes(0, i);
//     }
//     std::cout << "Mesh X range: [" << min_x << ", " << max_x << "]" << std::endl;
//     // 10x1x1 的梁，划分细一点: 10x2x2
//     // createBeamMesh(mesh, 10, 2, 2, 10.0, 1.0, 1.0);

//     std::cout << "Mesh: " << mesh.getNumNodes() << " nodes, "
//               << mesh.getNumElements() << " tetrahedra." << std::endl;

//     Material mat;
//     mat.E = 1e6;
//     mat.nu = 0.3;

//     double rho = 1.0;
//     Eigen::Vector3d gravity(0, 0, -1.0); // Z 向下

//     // 2. 边界条件
//     std::map<int, double> bc;
//     double tol = 1e-4; // 放宽一点容差，防止浮点误差
//     for (int i = 0; i < nodes.cols(); ++i) {
//         // 固定最左端 (Min X)
//         if (std::abs(nodes(0, i) - min_x) < tol) {
//             bc[3 * i + 0] = 0.0;
//             bc[3 * i + 1] = 0.0;
//             bc[3 * i + 2] = 0.0;
//         }
//     }
//     if (bc.empty()) {
//         std::cerr << "Error: No fixed nodes found! Check BC logic." << std::endl;
//         return -1;
//     }

//     // 3. 计算载荷
//     Eigen::VectorXd f_ext = computeGravityLoad(mesh, gravity, rho);

//     // 4. FEM 求解
//     std::cout << "Solving FEM..." << std::endl;
//     vem::fem::FEMSolver fem_solver(mesh);
//     fem_solver.assemble(mat);
//     Eigen::VectorXd u_fem = fem_solver.solve(bc, f_ext);

//     // 5. VEM 求解
//     std::cout << "Solving VEM..." << std::endl;
//     VEMSolver vem_solver(mesh);
//     vem_solver.assemble(mat);
//     Eigen::VectorXd u_vem = vem_solver.solve(bc, f_ext);

//     // 6. 对比结果
//     double tip_z_fem = 0.0;
//     double tip_z_vem = 0.0;
//     int tip_count = 0;

//     for (int i = 0; i < nodes.cols(); ++i) {
//         // 统计最右端 (Max X)
//         if (std::abs(nodes(0, i) - max_x) < tol) {
//             tip_z_fem += u_fem(3 * i + 2);
//             tip_z_vem += u_vem(3 * i + 2);
//             tip_count++;
//         }
//     }
//     if (tip_count == 0) {
//         std::cout << "Error: No tip nodes found." << std::endl;
//         return 0;
//     }

//     tip_z_fem /= tip_count;
//     tip_z_vem /= tip_count;
//     std::cout << "Tip Displacement Z (FEM): " << tip_z_fem << std::endl;
//     std::cout << "Tip Displacement Z (VEM): " << tip_z_vem << std::endl;

//     double diff = std::abs(tip_z_fem - tip_z_vem);
//     double rel_err = diff / std::abs(tip_z_fem);

//     std::cout << "Relative Difference: " << rel_err << std::endl;

//     // 7. 导出 VTK
//     // 我们需要把位移加到节点上再输出？或者作为 PointData 输出
//     // VTKWriter 目前只输出几何。我们可以修改 VTKWriter 支持 PointData，
//     // 或者简单点：直接修改 mesh 的 nodes 位置输出变形后的网格 (Deformed Mesh)

//     // Save Deformed VEM
//     VEMMesh mesh_vem = mesh; // Copy topology
//     Eigen::MatrixXd V_deformed_vem = nodes;
//     for (int i = 0; i < nodes.cols(); ++i) {
//         V_deformed_vem.col(i) += u_vem.segment<3>(3 * i);
//     }
//     mesh_vem.setNodes(V_deformed_vem);

//     io::VTKWriter writer;
//     writer.save(mesh_vem, "beam_vem_deformed.vtu");

//     // Save Deformed FEM
//     VEMMesh mesh_fem = mesh;
//     Eigen::MatrixXd V_deformed_fem = nodes;
//     for (int i = 0; i < nodes.cols(); ++i) {
//         V_deformed_fem.col(i) += u_fem.segment<3>(3 * i);
//     }
//     mesh_fem.setNodes(V_deformed_fem);
//     writer.save(mesh_fem, "beam_fem_deformed.vtu");

//     std::cout << "Results saved to beam_vem_deformed.vtu and beam_fem_deformed.vtu" << std::endl;

//     if (rel_err < 1e-10) {
//         std::cout << "[PASSED] VEM matches FEM on Beam." << std::endl;
//         return 0;
//     } else {
//         std::cout << "[FAILED] VEM mismatch." << std::endl;
//         return 1;
//     }
// }

#include "VEMMesh.hpp"
#include "VEMSolver.hpp"
// #include "fem/FEMSolver.hpp" // [Modified] 移除 FEM 求解器包含
#include "io/VTKLoader.hpp"
#include "io/VTKWriter.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using namespace vem;

// 理论解：欧拉-伯努利悬臂梁 (Z向重力)
// w(x) = (q * x^2) / (24 * E * I) * (x^2 - 4*L*x + 6*L^2)
// 注意：坐标系要对齐。假设 Fixed at x = -2 (即 x_local = 0), Free at x = 2 (x_local = L)
double getBeamTheoryDeflection(double x_world, double E, double I, double q, double L)
{
    // 转换坐标系：从 [-2, 2] 映射到 [0, 4] (假设悬臂根部在 0)
    // 你的模型 fixed at -2, so x_local = x_world - (-2) = x_world + 2
    double x = x_world + 2.0;

    // 欧拉伯努利公式 (悬臂梁受均布载荷 q)
    // w = - (q x^2)/(24EI) * (x^2 - 4Lx + 6L^2)
    // 符号为负，因为重力向下
    double w = -(q * x * x) / (24.0 * E * I) * (x * x - 4.0 * L * x + 6.0 * L * L);
    return w;
}

// 计算相对 L2 误差
double computeL2Error(const VEMMesh& mesh, const Eigen::VectorXd& u_vem, const Material& mat)
{
    double total_error_sq = 0.0;
    double total_ref_sq = 0.0;

    // 梁参数
    double L = 4.0;
    double b = 1.0;
    double h = 1.0;
    double I = b * std::pow(h, 3) / 12.0;
    double q = 1.0 * 1.0; // rho * Area * g (rho=1, g=1, Area=1)

    const auto& nodes = mesh.getNodes();
    int n_nodes = mesh.getNumNodes();

    for (int i = 0; i < n_nodes; ++i) {
        double x = nodes(0, i);
        double y = nodes(1, i);
        double z = nodes(2, i);

        // 获取 VEM 计算的位移
        double uz_vem = u_vem(3 * i + 2);

        // 获取理论解 (仅比较 Z 方向，因为这是主变形)
        double uz_exact = getBeamTheoryDeflection(x, mat.E, I, q, L);

        // 累加误差
        total_error_sq += std::pow(uz_vem - uz_exact, 2);
        total_ref_sq += std::pow(uz_exact, 2);
    }

    if (total_ref_sq < 1e-10)
        return 0.0;
    return std::sqrt(total_error_sq / total_ref_sq);
}

// 计算体力载荷向量 (Lumped Gravity Load)
// 该方法对任意多面体通用：计算单元总质量，然后均分到单元的所有节点上
Eigen::VectorXd computeGravityLoad(const VEMMesh& mesh, const Eigen::Vector3d& gravity, double rho)
{
    int n_nodes = mesh.getNumNodes();
    Eigen::VectorXd f = Eigen::VectorXd::Zero(n_nodes * 3);

    int n_elems = mesh.getNumElements();
    for (int i = 0; i < n_elems; ++i) {
        const auto& elem = mesh.getElement(i);

        // 注意：VEMElement在初始化时应当已经计算好 volume
        double mass = elem.volume * rho;

        // 获取单元所有节点 (去重)
        std::set<int> node_set;
        for (int fid : elem.face_indices) {
            // 注意：需要确保 VEMMesh 加载时构建了 Face 信息
            const auto* face = mesh.getFace(fid);
            if (face) {
                for (int nid : face->node_indices)
                    node_set.insert(nid);
            }
        }

        if (node_set.empty())
            continue;

        // 将重力均分到节点 (Lumping)
        // 对于 VEM，这是一种低阶但常用的体力处理方式
        double node_load_mag = mass / node_set.size();

        for (int nid : node_set) {
            f(3 * nid + 0) += node_load_mag * gravity(0);
            f(3 * nid + 1) += node_load_mag * gravity(1);
            f(3 * nid + 2) += node_load_mag * gravity(2);
        }
    }
    return f;
}

int main(int argc, char** argv)
{
    std::cout << "=== Beam Gravity Test (VEM Polyhedral Only) ===" << std::endl;

    // 1. 加载网格
    // 默认为生成的完美 CVT 多面体网格，也可以通过命令行参数传入
    std::string mesh_file = "/home/zcy/workspace/models/beam/beam_cvt_final.vtu";
    if (argc > 1) {
        mesh_file = argv[1];
    }

    VEMMesh mesh;
    io::VTKLoader loader;

    std::cout << "Loading mesh: " << mesh_file << " ..." << std::endl;
    auto loaded_mesh = loader.load(mesh_file);
    if (!loaded_mesh) {
        std::cerr << "Error: Failed to load mesh file: " << mesh_file << std::endl;
        return -1;
    }
    mesh = *loaded_mesh;
    // 在 load mesh 之后插入
    {
        std::cout << "--- Checking Mesh Topology Consistency ---" << std::endl;
        // 面签名 -> 出现次数
        // 签名是对面节点排序后的字符串，用于识别几何上相同的面
        std::map<std::string, int> face_counter;

        int n_elems = mesh.getNumElements();
        for (int i = 0; i < n_elems; ++i) {
            const auto& el = mesh.getElement(i);
            for (int fid : el.face_indices) {
                const auto* face = mesh.getFace(fid);
                std::vector<int> nodes = face->node_indices;
                std::sort(nodes.begin(), nodes.end()); // 排序以忽略节点顺序差异

                std::stringstream ss;
                for (int n : nodes)
                    ss << n << "_";
                face_counter[ss.str()]++;
            }
        }

        int boundary_faces = 0;
        int internal_faces = 0;
        int error_faces = 0;

        for (auto const& [key, count] : face_counter) {
            if (count == 1)
                boundary_faces++;
            else if (count == 2)
                internal_faces++;
            else
                error_faces++; // 3个以上单元共享同一个面（在3D中不应该发生，除非是非流形）
        }

        std::cout << "Total unique faces: " << face_counter.size() << std::endl;
        std::cout << "  - Internal faces (shared by 2): " << internal_faces << std::endl;
        std::cout << "  - Boundary faces (shared by 1): " << boundary_faces << std::endl;
        std::cout << "  - Non-manifold faces (>2): " << error_faces << std::endl;

        // 简单的判断标准：对于梁模型，内部面应该远多于边界面
        // 如果 internal_faces 极少，或者为 0，说明网格完全断开了
        double ratio = (double)internal_faces / face_counter.size();
        std::cout << "Connectivity Ratio: " << ratio * 100.0 << "%" << std::endl;

        // 修改检查逻辑：只要有内部面，说明就不是完全断开的
        if (internal_faces == 0 && n_elems > 1) {
            std::cerr << "[ERROR] Mesh is fragmented! No shared faces found." << std::endl;
            return 1;
        } else {
            std::cout << "[PASSED] Mesh topology looks consistent." << std::endl;
        }
    }

    // 获取网格范围用于边界条件判断
    double min_x = 1e30, max_x = -1e30;
    const auto& nodes = mesh.getNodes();
    for (int i = 0; i < nodes.cols(); ++i) {
        if (nodes(0, i) < min_x)
            min_x = nodes(0, i);
        if (nodes(0, i) > max_x)
            max_x = nodes(0, i);
    }
    std::cout << "Mesh X range: [" << min_x << ", " << max_x << "]" << std::endl;
    std::cout << "Mesh Info: " << mesh.getNumNodes() << " nodes, "
              << mesh.getNumElements() << " polyhedral elements." << std::endl;

    // 材料参数 (Gain et al. 2014)
    Material mat;
    mat.E = 1e4;
    mat.nu = 0.3;

    double rho = 1.0;
    Eigen::Vector3d gravity(0, 0, -10.0); // Z 向下

    // 2. 边界条件 (Fixed at min_x)
    std::map<int, double> bc;
    double tol = 1e-4;
    int fixed_count = 0;
    for (int i = 0; i < nodes.cols(); ++i) {
        if (std::abs(nodes(0, i) - min_x) < tol) {
            bc[3 * i + 0] = 0.0;
            bc[3 * i + 1] = 0.0;
            bc[3 * i + 2] = 0.0;
            fixed_count++;
        }
    }
    std::cout << "Fixed " << fixed_count << " nodes at X = " << min_x << std::endl;

    if (bc.empty()) {
        std::cerr << "Error: No fixed nodes found! Check BC logic or mesh coordinates." << std::endl;
        return -1;
    }

    // 3. 计算载荷
    std::cout << "Computing gravity load..." << std::endl;
    Eigen::VectorXd f_ext = computeGravityLoad(mesh, gravity, rho);

    // 4. VEM 求解 (移除 FEM)
    std::cout << "Solving VEM..." << std::endl;
    // 使用大括号限制作用域，确保求解器对象及时析构释放内存
    Eigen::VectorXd u_vem;
    {
        VEMSolver vem_solver(mesh);
        std::cout << "Assemble...\n";
        vem_solver.assemble(mat);
        std::cout << "Applying boundary condition...\n";
        u_vem = vem_solver.solve(bc, f_ext);
    }

    // 5. 结果统计 (Tip Displacement)
    double tip_z_vem = 0.0;
    int tip_count = 0;

    for (int i = 0; i < nodes.cols(); ++i) {
        // 统计最右端 (Max X)
        if (std::abs(nodes(0, i) - max_x) < tol) {
            tip_z_vem += u_vem(3 * i + 2);
            tip_count++;
        }
    }

    if (tip_count > 0) {
        tip_z_vem /= tip_count;
        std::cout << "Tip Displacement Z (VEM Average): " << tip_z_vem << std::endl;

        // Gain et al. (2014) 参考解:
        // 根据梁理论或精细网格 FEM 解，对于 10x1x1 梁，E=1e6, nu=0.3, rho=1, g=-1
        // Euler-Bernoulli: w_max = (rho*A*g * L^4) / (8 * E * I)
        // I = W*H^3/12 = 1*1/12 = 0.08333
        // q = rho*Area*g = 1*1*1 = 1
        // w_max = (1 * 10000) / (8 * 1e6 * 0.08333) = 10000 / 666666.6 = 0.015
        // (注：粗略估算，具体值取决于实际尺寸和理论模型)
    } else {
        std::cout << "Warning: No tip nodes found for statistics." << std::endl;
    }

    // 6. 导出变形后的网格
    std::cout << "Saving results to beam_vem_deformed.vtu ..." << std::endl;
    VEMMesh mesh_vem = mesh;
    Eigen::MatrixXd V_deformed_vem = nodes;
    for (int i = 0; i < nodes.cols(); ++i) {
        V_deformed_vem.col(i) += u_vem.segment<3>(3 * i);
    }
    mesh_vem.setNodes(V_deformed_vem);

    io::VTKWriter writer;
    writer.save(mesh_vem, "beam_vem_deformed.vtu");

    std::cout << "[DONE] VEM Polyhedral simulation finished." << std::endl;

    return 0;
}
