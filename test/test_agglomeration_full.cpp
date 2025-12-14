#include "AgglomerationManager.hpp"
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"
#include "io/VTKLoader.hpp"
#include "io/VTKWriter.hpp"
#include <iostream>
#include <map> // 必须包含
#include <set>

using namespace vem; // 使用命名空间简化代码

// 计算重力载荷
void computeGravityLoad(const VEMMesh& mesh, Eigen::VectorXd& f)
{
    int n_nodes = mesh.getNumNodes();
    f.setZero(3 * n_nodes);

    Eigen::Vector3d g(0, 0, -1.0);
    double rho = 1.0;

    for (int i = 0; i < mesh.getNumElements(); ++i) {
        if (!mesh.isElementActive(i))
            continue;
        const auto& el = mesh.getElement(i);

        double mass = el.volume * rho;

        std::set<int> unique_nodes;
        for (int fid : el.face_indices) {
            const auto* face = mesh.getFace(fid);
            for (int nid : face->node_indices)
                unique_nodes.insert(nid);
        }

        if (unique_nodes.empty())
            continue;

        Eigen::Vector3d force_per_node = (mass * g) / double(unique_nodes.size());
        for (int nid : unique_nodes) {
            f.segment<3>(3 * nid) += force_per_node;
        }
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./test_agglomeration_full <mesh_file>" << std::endl;
        return 1;
    }

    // 1. 加载网格
    std::string mesh_file = argv[1];
    io::VTKLoader loader;
    auto mesh_ptr = loader.load(mesh_file);
    if (!mesh_ptr)
        return 1;
    VEMMesh& mesh = *mesh_ptr;

    const auto& nodes = mesh.getNodes();
    double min_x = nodes.row(0).minCoeff();
    std::cout << "Mesh Bounds X min: " << min_x << std::endl;

    Material mat;
    mat.E = 1e6;
    mat.nu = 0.3;

    // 2. 执行 Agglomeration
    AgglomerationManager agglomerator(mesh, mat);
    agglomerator.run(0.1, 5);

    // 3. 垃圾回收
    std::cout << "Cleaning up mesh..." << std::endl;
    mesh.garbageCollect();

    // 4. 求解
    std::cout << "Solving on optimized mesh..." << std::endl;
    VEMSolver solver(mesh);

    // [FIX] 组装刚度矩阵
    solver.assemble(mat);

    // [FIX] 构建 Dirichlet 边界条件 Map (Global DOF Index -> Value)
    std::map<int, double> dirichlet_bc;
    double tol = 1e-3;
    int fixed_count = 0;

    for (int i = 0; i < mesh.getNumNodes(); ++i) {
        if (std::abs(nodes(0, i) - min_x) < tol) {
            // 固定该节点的所有 3 个自由度 (x, y, z)
            dirichlet_bc[3 * i + 0] = 0.0;
            dirichlet_bc[3 * i + 1] = 0.0;
            dirichlet_bc[3 * i + 2] = 0.0;
            fixed_count++;
        }
    }
    std::cout << "Fixed " << fixed_count << " nodes at X = " << min_x << std::endl;

    // Load
    Eigen::VectorXd f;
    computeGravityLoad(mesh, f);

    // Solve
    Eigen::VectorXd u = solver.solve(dirichlet_bc, f);

    // 5. 保存结果
    io::VTKWriter writer;
    std::string out_name = "beam_agglomerated_result.vtu";
    writer.save(mesh, u, out_name);
    std::cout << "Result saved to " << out_name << std::endl;

    return 0;
}
