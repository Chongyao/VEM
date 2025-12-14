#include "AgglomerationManager.hpp"
#include "VEMMesh.hpp"
#include "io/VTKLoader.hpp"
#include <iomanip>
#include <iostream>

using namespace vem;

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./test_agglomeration_preview <mesh_file>" << std::endl;
        return 1;
    }

    std::string mesh_file = argv[1];
    io::VTKLoader loader;
    auto mesh_ptr = loader.load(mesh_file);
    if (!mesh_ptr)
        return 1;
    VEMMesh& mesh = *mesh_ptr;

    std::cout << "Mesh loaded: " << mesh.getNumElements() << " elements." << std::endl;

    Material mat;
    mat.E = 1e6;
    mat.nu = 0.3;
    AgglomerationManager algo(mesh, mat);

    // 1. 找出当前最差的一个单元 (复用之前的逻辑或简单遍历)
    int worst_id = -1;
    double min_sigma = 1.0;

    std::cout << "Finding worst element..." << std::endl;
    for (int i = 0; i < mesh.getNumElements(); ++i) {
        double s = algo.computeCurrentSigma(i);
        if (s < min_sigma) {
            min_sigma = s;
            worst_id = i;
        }
    }

    std::cout << "\n=== Worst Element Analysis ===" << std::endl;
    std::cout << "Worst Element ID: " << worst_id << std::endl;
    std::cout << "Current Sigma:    " << std::scientific << min_sigma << std::endl;
    std::cout << "Neighbors:        " << mesh.getElementNeighbors(worst_id).size() << std::endl;

    std::cout << "\n--- Simulation: Merging with Neighbors ---" << std::endl;
    std::cout << std::left << std::setw(15) << "Neighbor ID"
              << std::setw(15) << "New Sigma"
              << std::setw(15) << "Improvement" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    std::vector<int> neighbors = mesh.getElementNeighbors(worst_id);
    int best_nb = -1;
    double best_s = -1.0;

    for (int nb : neighbors) {
        double new_s = algo.simulateMergeSigma(worst_id, nb);
        double imp = new_s / min_sigma;

        std::cout << std::left << std::setw(15) << nb
                  << std::setw(15) << std::scientific << std::setprecision(4) << new_s
                  << std::setw(15) << std::fixed << std::setprecision(2) << imp << "x" << std::endl;

        if (new_s > best_s) {
            best_s = new_s;
            best_nb = nb;
        }
    }

    std::cout << "\n[Recommendation] Merge Element " << worst_id << " with " << best_nb
              << " -> New Sigma: " << best_s
              << " (" << (best_s / min_sigma) << "x improvement)" << std::endl;

    return 0;
}
