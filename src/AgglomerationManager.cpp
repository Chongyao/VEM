#include "AgglomerationManager.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>

// 辅助：获取合并后的面列表 (对称差集)
std::vector<int> AgglomerationManager::getMergedFaces(int id1, int id2)
{
    const auto& el1 = mesh_.getElement(id1);
    const auto& el2 = mesh_.getElement(id2);

    std::set<int> f1(el1.face_indices.begin(), el1.face_indices.end());
    std::set<int> f2(el2.face_indices.begin(), el2.face_indices.end());

    std::vector<int> new_faces;
    // 保留不共享的面
    for (int f : f1) {
        if (f2.find(f) == f2.end()) {
            new_faces.push_back(f);
        }
    }
    for (int f : f2) {
        if (f1.find(f) == f1.end()) {
            new_faces.push_back(f);
        }
    }

    return new_faces;
}

double AgglomerationManager::computeCurrentSigma(int elem_id)
{
    if (!mesh_.isElementActive(elem_id)) {
        return 0.0;
    }

    const auto& elem = mesh_.getElement(elem_id);
    VEMElement vem_el(mesh_, elem);
    Eigen::MatrixXd K = vem_el.computeStiffness(mat_);

    return MathUtils::computeStabilityRatio(K);
}

double AgglomerationManager::simulateMergeSigma(int id1, int id2)
{
    // 1. 获取新面集合
    std::vector<int> new_faces = getMergedFaces(id1, id2);
    if (new_faces.empty()) {
        return 0.0;
    }

    // 2. 构造临时单元
    PolyhedronElement temp_el;
    temp_el.id = -1;
    temp_el.face_indices = new_faces;
    temp_el.active = true;

    // 3. 准备几何计算所需指针
    std::vector<const VEMFace*> face_ptrs;
    face_ptrs.reserve(new_faces.size());
    for (int fid : new_faces) {
        face_ptrs.push_back(mesh_.getFace(fid));
    }

    // 4. 调用 GeometryUtils 计算属性 (复用核心逻辑)
    GeometryUtils::computePolyhedronProps(mesh_.getNodes(), face_ptrs, temp_el.volume, temp_el.centroid);

    // 5. 计算刚度与 Sigma
    if (temp_el.volume <= 1e-12) {
        return 0.0; // 避免退化单元
    }

    VEMElement vem_el(mesh_, temp_el);
    Eigen::MatrixXd K = vem_el.computeStiffness(mat_);

    return MathUtils::computeStabilityRatio(K);
}

int AgglomerationManager::mergeElements(int id1, int id2)
{
    // 1. 计算新面
    std::vector<int> new_faces = getMergedFaces(id1, id2);

    // 2. 添加新单元 (Mesh 会自动计算几何、更新拓扑)
    int new_id = mesh_.addElement(new_faces);

    // 3. 标记旧单元无效
    mesh_.deactivateElement(id1);
    mesh_.deactivateElement(id2);

    return new_id;
}

MergeCandidate AgglomerationManager::findBestNeighbor(int elem_id)
{
    std::vector<int> neighbors = mesh_.getElementNeighbors(elem_id);

    int best_id = -1;
    double max_sigma = -1.0;

    // 当前自己的 sigma，用于判断合并是否有提升
    // double current_sigma = computeCurrentSigma(elem_id);

    for (int nb_id : neighbors) {
        if (!mesh_.isElementActive(nb_id)) {
            continue;
        }

        double s = simulateMergeSigma(elem_id, nb_id);

        // 简单的贪心策略：找合并后 sigma 最大的
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

void AgglomerationManager::run(double sigma_threshold, int max_passes)
{
    std::cout << "\n=== Starting Agglomeration Process ===\n";
    std::cout << "Threshold: " << sigma_threshold << ", Max Passes: " << max_passes << "\n";

    for (int pass = 0; pass < max_passes; ++pass) {
        int merge_count = 0;
        int active_count = 0;

        // 1. 构建优先队列 (只包含活跃的坏单元)
        // 使用 std::priority_queue 需要 vector 容器和 greater 比较器
        std::priority_queue<ElementEntry, std::vector<ElementEntry>, std::greater<ElementEntry>> pq;

        for (int i = 0; i < mesh_.getNumElements(); ++i) {
            if (!mesh_.isElementActive(i)) {
                continue;
            }
            active_count++;

            double sigma = computeCurrentSigma(i);
            if (sigma < sigma_threshold) {
                pq.push({ i, sigma });
            }
        }

        std::cout << "Pass " << pass + 1 << ": Active Elements = " << active_count
                  << ", Bad Elements = " << pq.size() << "\n";

        if (pq.empty()) {
            std::cout << "-> No bad elements left. Converged." << "\n";
            break;
        }

        // 2. 处理队列
        // 注意：队列中的单元可能在处理前面的单元时被合并了(作为邻居)，所以要再次检查 active
        while (!pq.empty()) {
            ElementEntry entry = pq.top();
            pq.pop();
            int curr_id = entry.id;

            // 检查是否仍然活跃
            if (!mesh_.isElementActive(curr_id)) {
                continue;
            }

            // 寻找合并对象
            MergeCandidate best = findBestNeighbor(curr_id);

            // 接受合并的条件 (CutVEM paper criterion)
            // sigma_new > sigma_old (或者稍微放宽)
            // 这里我们严格一点：必须有所提升

            if (best.neighbor_id != -1) {
                // 检查提升幅度
                // 论文建议: sigma_new > min(threshold, beta * sigma_old)
                // 这里简单处理：只要比原来好，或者达到了阈值就可以

                // 如果合并后依然很差且没有提升，可能就不合并了？
                // 我们的策略：只要 max_sigma > entry.sigma * 1.01 (1%提升) 就合并
                // 或者 max_sigma > threshold

                if (best.new_sigma > entry.sigma || best.new_sigma > sigma_threshold) {
                    // 执行合并
                    int new_id = mergeElements(curr_id, best.neighbor_id);
                    merge_count++;

                    // 新单元暂不加入本轮队列 (Wait for next pass)
                    // 这是为了防止单次 pass 产生超级大单元
                }
            }
        }

        std::cout << "-> Merged " << merge_count << " pairs." << std::endl;
        if (merge_count == 0) {
            std::cout << "-> Stalled (could not find better neighbors). Stopping." << std::endl;
            break;
        }
    }

    // 最终统计
    int final_active = 0;
    for (int i = 0; i < mesh_.getNumElements(); ++i) {
        if (mesh_.isElementActive(i)) {
            final_active++;
        }
    }

    std::cout << "=== Agglomeration Finished. Final Active Elements: " << final_active << " ===\n";
}
