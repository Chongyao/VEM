#pragma once
#include "MathUtils.hpp"
#include "VEMElement.hpp"
#include "VEMMesh.hpp"
#include <map>
#include <set>
#include <vector>

namespace vem {

// 简单的合并候选项结构
struct MergeCandidate {
    int neighbor_id;
    double new_sigma;
};

class AgglomerationManager {
private:
    VEMMesh& mesh_;
    Material mat_;

    // 标记单元是否已被合并 (active = true 表示有效)
    std::vector<bool> active_mask_;

public:
    AgglomerationManager(VEMMesh& mesh, Material mat)
        : mesh_(mesh)
        , mat_(mat)
    {
        // 初始化掩码，默认全有效
        active_mask_.resize(mesh_.getNumElements(), true);
    }

    // Step 1 核心功能：模拟合并
    // 计算如果将 el1 和 el2 合并，新单元的 sigma 是多少
    // 返回: 0.0 ~ 1.0 的比率
    double simulateMergeSigma(int id1, int id2);

    // Step 1 核心功能：寻找最佳邻居
    // 遍历所有邻居，找到合并后 sigma 最高的那个
    // 返回: {best_neighbor_id, max_sigma}。如果没有合法合并，返回 {-1, 0.0}
    MergeCandidate findBestNeighbor(int elem_id);

    // 计算当前单元的 sigma (工具函数)
    double computeCurrentSigma(int elem_id);

private:
    // 辅助：为临时单元计算几何属性 (Volume, Centroid)
    // 因为 PolyhedronElement::updateGeometricProps 可能需要访问私有 faces 数组，
    // 我们在这里实现一个适配版本
    void updateTempElementGeometry(PolyhedronElement& elem);
};

} // namespace vem
