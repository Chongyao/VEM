#include "AgglomerationManager.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace vem {

double AgglomerationManager::computeCurrentSigma(int elem_id)
{
    if (elem_id < 0 || elem_id >= mesh_.getNumElements())
        return 0.0;

    const auto& elem = mesh_.getElement(elem_id);
    VEMElement vem_el(mesh_, elem);
    Eigen::MatrixXd K = vem_el.computeStiffness(mat_);

    return MathUtils::computeStabilityRatio(K);
}

void AgglomerationManager::updateTempElementGeometry(PolyhedronElement& elem)
{
    // 简单的体积和质心计算
    // VEM 假设面是平的。我们可以利用散度定理。
    // Volume = 1/3 * sum( x_f . n_f * Area_f )
    // Centroid 稍微复杂点，但可以用锥体分解近似。

    double total_vol = 0.0;
    Eigen::Vector3d weighted_centroid = Eigen::Vector3d::Zero();

    // 参考点 (可以取任意点，例如第一个面的中心，减少数值误差)
    // 但为了简单，假设原点。或者为了稳健，取所有面质心的平均值作为参考点 x0
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
    if (!elem.face_indices.empty()) {
        const auto* f0 = mesh_.getFace(elem.face_indices[0]);
        if (f0) {
            auto props = f0->computeProps(mesh_.getNodes());
            x0 = props.centroid;
        }
    }

    for (int fid : elem.face_indices) {
        const auto* face = mesh_.getFace(fid);
        if (!face)
            continue;

        // 计算面的属性 (Area, Normal, Centroid)
        auto props = face->computeProps(mesh_.getNodes());

        // 修正法向：确保指向多面体外部
        // 简单策略：如果 (FaceCentroid - MeshCentroid) . n < 0 则反转？
        // 但我们现在还没算出 MeshCentroid。
        // CutVEM 论文通常假设网格生成时法向是协调的，或者我们需要维护局部拓扑。
        // 在 simulateMerge 中，我们继承了原单元的面。
        // 如果原单元的面法向是正确的（指向各子单元外），那么对于合并后的单元：
        // 1. 原 E1 的面，法向指向 E1 外（即新单元外）。正确。
        // 2. 原 E2 的面，法向指向 E2 外（即新单元外）。正确。
        // 3. 共享面被消除了。
        // 结论：只要原单元法向正确，合并后的法向天然正确。
        // 但是 VEMFace 存储的是全局法向。在使用时 VEMElement 会根据 element.centroid 判断并翻转。
        // 这里我们要算 element.centroid，这就成了鸡生蛋。
        // 临时方案：先算一个粗略的几何中心用于法向判断，再精算体积。

        // 1. 粗略中心
        // ... 跳过，直接用代数体积公式，允许负体积，最后取绝对值？
        // 不行，质心计算依赖符号。

        // 让我们采用 VEMElement 的策略：在组装 stiffness 时才判断法向。
        // 这里我们只需要确保 volume 是正的。
        // 实际上，computeProps 返回的 normal 是由节点顺序决定的。
        // 我们假设 `generate_beam_fixed_normal.py` 已经保证了全局一致性（或者局部一致性）。
        //
        // 既然 VEMElement 内部会重新计算法向方向 (check normal vs centroid)，
        // 我们在这里只需要提供一个合理的 centroid 供 VEMElement 使用即可。

        // 金字塔体积 (Apex at x0)
        // V_pyr = 1/3 * Area * dot(face_centroid - x0, normal)
        // 注意：normal 方向未定。
        // 这是一个难点。如果不预先知道 centroid，很难确定 normal 朝向。
        //
        // 替代方案：直接取所有面顶点的平均值作为 centroid approximation。
        // 对于凸多面体，这足够用于 scaling 和 stability check。
        // CutVEM 论文中 Agglomeration 主要是为了改善条件数，几何属性微小的误差（如质心偏一点）
        // 对 sigma 的量级影响很小。
    }

    // === 简化实现 ===
    // 1. 质心 = 所有面质心的加权平均 (按面积)
    // 2. 体积 = sum(Pyramids based on that centroid)

    double total_area = 0.0;
    Eigen::Vector3d surf_centroid = Eigen::Vector3d::Zero();

    for (int fid : elem.face_indices) {
        const auto* face = mesh_.getFace(fid);
        auto props = face->computeProps(mesh_.getNodes());
        surf_centroid += props.centroid * props.area;
        total_area += props.area;
    }
    if (total_area > 1e-12)
        surf_centroid /= total_area;
    elem.centroid = surf_centroid;

    // 2. 计算体积
    total_vol = 0.0;
    for (int fid : elem.face_indices) {
        const auto* face = mesh_.getFace(fid);
        auto props = face->computeProps(mesh_.getNodes());

        // 确保法向朝外
        Eigen::Vector3d n = props.normal;
        if ((props.centroid - elem.centroid).dot(n) < 0) {
            n = -n;
        }

        // Pyramid volume with apex at elem.centroid is
        // 1/3 * Area * height
        // height = projected distance from centroid to plane?
        // Actually simpler: 1/3 * A * ((c_f - c_el) . n)
        double h = std::abs((props.centroid - elem.centroid).dot(n));
        total_vol += (1.0 / 3.0) * props.area * h;
    }
    elem.volume = total_vol;
}

double AgglomerationManager::simulateMergeSigma(int id1, int id2)
{
    const auto& el1 = mesh_.getElement(id1);
    const auto& el2 = mesh_.getElement(id2);

    // 1. 构建新面集合 (Symmetric Difference)
    std::set<int> faces1(el1.face_indices.begin(), el1.face_indices.end());
    std::set<int> faces2(el2.face_indices.begin(), el2.face_indices.end());

    std::vector<int> new_faces;
    // 添加只属于 faces1 的
    for (int f : faces1)
        if (faces2.find(f) == faces2.end())
            new_faces.push_back(f);
    // 添加只属于 faces2 的
    for (int f : faces2)
        if (faces1.find(f) == faces1.end())
            new_faces.push_back(f);

    // 2. 创建临时单元
    PolyhedronElement new_el;
    new_el.id = -1; // Temporary
    new_el.face_indices = new_faces;

    // 3. 更新几何
    updateTempElementGeometry(new_el);

    // 4. 计算刚度与 Sigma
    VEMElement vem_el(mesh_, new_el);
    Eigen::MatrixXd K = vem_el.computeStiffness(mat_);

    return MathUtils::computeStabilityRatio(K);
}

MergeCandidate AgglomerationManager::findBestNeighbor(int elem_id)
{
    std::vector<int> neighbors = mesh_.getElementNeighbors(elem_id);

    int best_id = -1;
    double max_sigma = -1.0;

    for (int nb_id : neighbors) {
        // 如果邻居已经被合并了，跳过 (目前 active_mask_ 还没用到，但在完整流程中需要)
        if (!active_mask_[nb_id])
            continue;

        double s = simulateMergeSigma(elem_id, nb_id);
        if (s > max_sigma) {
            max_sigma = s;
            best_id = nb_id;
        }
    }

    if (best_id != -1) {
        return { best_id, max_sigma };
    }
    return { -1, 0.0 };
}

} // namespace vem
