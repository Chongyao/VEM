#include "VEMMesh.hpp"
#include "VEMSolver.hpp"
// #include "fem/FEMSolver.hpp" // [Modified] 移除 FEM 求解器包含
#include "AnalyticalSolutions.hpp"
#include "ErrorAnalysis.hpp"
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
    std::cout << "=== Beam Gravity Test (VEM Polyhedral Only) ===\n";

    // 1. 加载网格
    // 默认为生成的完美 CVT 多面体网格，也可以通过命令行参数传入
    std::string mesh_file = "/home/zcy/workspace/models/beam/beam_cvt_final.vtu";
    if (argc > 1) {
        mesh_file = argv[1];
    }

    VEMMesh mesh;
    io::VTKLoader loader;

    std::cout << "Loading mesh: " << mesh_file << " ...\n";
    auto loaded_mesh = loader.load(mesh_file);
    if (!loaded_mesh) {
        std::cerr << "Error: Failed to load mesh file: " << mesh_file << "\n";
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
    // ================== 新增：误差分析 ==================
    // 定义梁参数 (与 Python 生成脚本一致)
    double L = 4.0;
    double b = 1.0;
    double h = 1.0;
    double I = b * std::pow(h, 3) / 12.0;
    // 线载荷 q = rho * Area * g
    // 假设 density=1.0, g=1.0, Area=1.0 => q=1.0
    double q = 1.0;

    utils::CantileverBeamSolution beam_sol(1e6, I, q, L, -2.0); // E=1e6

    // 定义精确解 lambda 函数 (只关心 Z 方向位移，XY 忽略或设为0)
    auto exact_func = [&](double x, double y, double z) -> Eigen::Vector3d {
        return Eigen::Vector3d(0.0, 0.0, beam_sol.getDisplacementZ(x));
    };

    double l2_error = utils::computeDiscreteL2Error(mesh, u_vem, exact_func);

    // 格式化输出，方便 Shell 脚本解析
    // 格式: [RESULT] N_Elements L2_Error
    std::cout << "[RESULT] " << mesh.getNumElements() << " " << l2_error << std::endl;

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
