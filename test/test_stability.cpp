#include "../src/MathUtils.hpp" // 引入新工具
#include "VEMElement.hpp"
#include "VEMMesh.hpp"
#include "io/VTKLoader.hpp"
#include <iomanip> // 格式化输出
#include <iostream>
#include <queue> // 优先队列
#include <vector>

using namespace vem;

// 定义一个简单的结构体来存储单元质量信息
struct ElementQuality {
    int id;
    double sigma;
    int n_faces;

    // 优先队列默认是最大堆（Top 是最大的）。
    // 我们希望 Top 是 sigma 最小的（最差的）。
    // 所以我们需要 "大于" 运算符：如果 a.sigma > b.sigma，则 a "大于" b。
    // std::priority_queue<T, vector<T>, greater<T>> 会让最小的元素浮在上面。
    // 这里我们重载 > 运算符来实现同样的效果。
    bool operator>(const ElementQuality& other) const
    {
        return sigma > other.sigma;
    }
};

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: ./test_stability <mesh_file>" << std::endl;
        return 1;
    }

    // 1. 加载网格
    std::string mesh_file = argv[1];
    io::VTKLoader loader;
    auto mesh_ptr = loader.load(mesh_file);
    if (!mesh_ptr) {
        std::cerr << "Failed to load mesh: " << mesh_file << std::endl;
        return 1;
    }
    VEMMesh& mesh = *mesh_ptr;

    std::cout << "Mesh loaded: " << mesh.getNumElements() << " elements." << std::endl;

    // 2. 准备材料 (仅用于刚度计算，不影响相对比率)
    Material mat;
    mat.E = 1e6;
    mat.nu = 0.3;

    // 3. 定义优先队列 (最小堆，sigma 小的在 top)
    // 使用 std::greater 和我们在结构体中定义的 operator> 配合?
    // 其实直接用 std::priority_queue<ElementQuality, std::vector<ElementQuality>, std::greater<ElementQuality>>
    // 需要重载 > 运算符。
    std::priority_queue<ElementQuality, std::vector<ElementQuality>, std::greater<ElementQuality>> pq;

    std::cout << "--- Computing Stability Ratios ---" << std::endl;

    for (int i = 0; i < mesh.getNumElements(); ++i) {
        const auto& elem = mesh.getElement(i);

        // A. 使用 VEMElement 计算刚度矩阵 (复用现有代码)
        VEMElement vem_el(mesh, elem);
        Eigen::MatrixXd K = vem_el.computeStiffness(mat);

        // B. 使用 MathUtils 计算稳定性比率
        double sigma = MathUtils::computeStabilityRatio(K);

        // C. 入队
        pq.push({ i, sigma, (int)elem.face_indices.size() });
    }

    // 4. 输出最差的 10 个单元 (模拟 CutVEM 的处理顺序)
    std::cout << "\n=== Top 10 Worst Elements (Candidates for Agglomeration) ===" << std::endl;
    std::cout << std::left << std::setw(10) << "Rank"
              << std::setw(10) << "ElemID"
              << std::setw(15) << "Sigma"
              << std::setw(10) << "Faces" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    int count = 0;
    while (!pq.empty() && count < 10) {
        ElementQuality bad = pq.top();
        pq.pop();

        std::cout << std::left << std::setw(10) << count + 1
                  << std::setw(10) << bad.id
                  << std::setw(15) << std::scientific << std::setprecision(4) << bad.sigma
                  << std::setw(10) << bad.n_faces << std::endl;
        count++;
    }

    return 0;
}
