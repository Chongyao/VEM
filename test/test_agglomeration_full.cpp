#include "AgglomerationManager.hpp"
#include "VEMMesh.hpp"
#include "VEMSolver.hpp"
#include "io/VTKLoader.hpp"
#include "io/VTKWriter.hpp"
#include <iostream>
#include <set>

using namespace vem;
// 计算重力载荷 (复用逻辑)
void computeGravityLoad(const VEMMesh& mesh, Eigen::VectorXd& f)
{
    int n_nodes = mesh.getNumNodes();
    f.setZero(3 * n_nodes);

    Eigen::Vector3d g(0, 0, -1.0); // 重力向下
    double rho = 1.0;

    for (int i = 0; i < mesh.getNumElements(); ++i) {
        const auto& el = mesh.getElement(i);
        // 注意：GC 之后所有剩下的单元都是 active 的，不需要检查 mask

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

    // 2. 获取几何范围 (Bounding Box)
    const auto& nodes = mesh.getNodes();
    double min_x = nodes.row(0).minCoeff();
    double max_x = nodes.row(0).maxCoeff();
    std::cout << "Mesh Bounds X: [" << min_x << ", " << max_x << "]" << std::endl;

    Material mat;
    mat.E = 1e6;
    mat.nu = 0.3;

    // 3. 执行 Agglomeration
    AgglomerationManager agglomerator(mesh, mat);

    // 设定阈值 (例如 0.1) 和最大轮数
    // 如果网格本身质量较好，可以把阈值设高一点来强制测试合并效果
    agglomerator.run(0.1, 5);

    // 4. [关键] 垃圾回收
    // 清除所有被合并的死单元，重排 ID，重建拓扑
    std::cout << "Cleaning up mesh..." << std::endl;
    mesh.garbageCollect();

    // 5. 求解 (在干净的网格上)
    std::cout << "Solving on optimized mesh (" << mesh.getNumElements() << " elements)..." << std::endl;

    VEMSolver solver(mesh);

    // BC: Fix X = min_x (动态范围)
    std::vector<int> fixed_nodes;
    double tol = 1e-3;
    for (int i = 0; i < mesh.getNumNodes(); ++i) {
        if (std::abs(nodes(0, i) - min_x) < tol) {
            fixed_nodes.push_back(i);
        }
    }
    std::cout << "Fixed " << fixed_nodes.size() << " nodes at X = " << min_x << std::endl;
    solver.setBoundaryConditions(fixed_nodes);

    // Load
    Eigen::VectorXd f;
    computeGravityLoad(mesh, f);

    // Solve
    Eigen::VectorXd u = solver.solve(f, mat);

    // 6. 保存结果
    // 由于执行了 GC，Writer 不需要做任何修改，直接遍历输出即可
    io::VTKWriter writer;
    std::string out_name = "beam_agglomerated_result.vtu";
    writer.save(mesh, u, out_name);
    std::cout << "Result saved to " << out_name << std::endl;

    return 0;
}
