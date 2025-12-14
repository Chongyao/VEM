#pragma once
#include "MathUtils.hpp"
#include "VEMElement.hpp"
#include "VEMMesh.hpp"
#include "VEMSolver.hpp" // 用于计算刚度
#include <queue>
#include <vector>

// 定义合并候选项
struct MergeCandidate {
    int neighbor_id;
    double new_sigma;
};

// 定义优先队列元素
struct ElementEntry {
    int id;
    double sigma;
    // 最小堆：sigma 小的排在前面
    bool operator>(const ElementEntry& other) const
    {
        return sigma > other.sigma;
    }
};

class AgglomerationManager {
private:
    VEMMesh& mesh_;
    Material mat_;

    // 我们直接使用 mesh_.isElementActive()，不再需要单独维护 mask
    // 但为了算法效率，可以用局部变量缓存

public:
    AgglomerationManager(VEMMesh& mesh, Material mat)
        : mesh_(mesh)
        , mat_(mat)
    {
    }

    // [New] 核心算法主循环
    // sigma_threshold: 低于此值的单元尝试合并 (例如 0.1)
    // max_passes: 最大迭代轮数
    void run(double sigma_threshold, int max_passes = 5);

    // [New] 执行合并：将 id1 和 id2 合并，返回新单元 ID
    int mergeElements(int id1, int id2);

    // 模拟合并：计算如果合并后的 sigma
    double simulateMergeSigma(int id1, int id2);

    // 寻找最佳邻居
    MergeCandidate findBestNeighbor(int elem_id);

    // 工具：计算当前单元 Sigma
    double computeCurrentSigma(int elem_id);

private:
    // 辅助：获取两个单元合并后的面集合
    std::vector<int> getMergedFaces(int id1, int id2);
};
