#include "VEMMesh.hpp"
#include "io/VTKLoader.hpp"
#include <iostream>
#include <vector>

using namespace vem;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./test_adjacency <mesh_file>" << std::endl;
        return 1;
    }

    std::string mesh_file = argv[1];
    io::VTKLoader loader;
    auto mesh_ptr = loader.load(mesh_file);
    if (!mesh_ptr)
        return 1;
    VEMMesh& mesh = *mesh_ptr;

    std::cout << "\n=== Testing Adjacency ===" << std::endl;

    // 1. 验证流形性
    if (!mesh.checkManifold()) {
        std::cerr << "[FAILED] Mesh is not 2-manifold (some faces shared by >2 elements)." << std::endl;
        return 1;
    }

    // 2. 随机抽查几个单元的邻居
    int n_elems = mesh.getNumElements();
    int samples = std::min(5, n_elems);

    for (int i = 0; i < samples; ++i) {
        // 比如取中间的几个 ID
        int id = (n_elems / (samples + 1)) * (i + 1);

        std::vector<int> neighbors = mesh.getElementNeighbors(id);

        std::cout << "Element " << id << " has " << neighbors.size() << " neighbors: ";
        for (int n : neighbors)
            std::cout << n << " ";
        std::cout << std::endl;

        // 简单一致性检查：如果 A 是 B 的邻居，B 必须是 A 的邻居
        for (int n : neighbors) {
            std::vector<int> neighbors_of_n = mesh.getElementNeighbors(n);
            bool found_back = false;
            for (int nn : neighbors_of_n)
                if (nn == id)
                    found_back = true;

            if (!found_back) {
                std::cerr << "[Error] Adjacency asymmetry: " << id << " sees " << n << ", but " << n << " does not see " << id << std::endl;
                return 1;
            }
        }
    }

    std::cout << "[PASSED] Topology graph is consistent." << std::endl;

    return 0;
}
