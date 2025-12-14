#include "VEMMesh.hpp"
#include <fstream>
#include <iostream>
VEMMesh::VEMMesh(const VEMMesh& other)
    : nodes_(other.nodes_)
    , elements_(other.elements_)
{
    faces_.reserve(other.faces_.size());
    for (const auto& face_ptr : other.faces_) {
        faces_.push_back(face_ptr->clone());
    }
}

VEMMesh& VEMMesh::operator=(const VEMMesh& other)
{
    if (this != &other) {
        nodes_ = other.nodes_;
        elements_ = other.elements_;

        faces_.clear();
        faces_.reserve(other.faces_.size());
        for (const auto& face_ptr : other.faces_) {
            faces_.push_back(face_ptr->clone());
        }
    }
    return *this;
}
void PolyhedronElement::updateGeometricProps(const Eigen::MatrixXd& nodes,
    const std::vector<std::unique_ptr<VEMFace>>& faces)
{
    volume = 0.0;
    Eigen::Vector3d moment_1st = Eigen::Vector3d::Zero();

    // 1. 计算参考点 X_ref (简单的顶点平均)
    // 这是一个鲁棒性技巧：对于凸多面体，任何点都可以；
    // 对于星形多面体，内核点最好；均值点通常落在内核中。
    Eigen::Vector3d x_ref = Eigen::Vector3d::Zero();
    int node_count = 0;
    // 这里为了简单，我们遍历面的所有点（会有重复，但对均值影响不大，或者可以用 std::set 去重）
    // 为了极致性能，后面可以优化，现在先求面的重心的平均
    for (int fid : face_indices) {
        // 使用面的质心来计算体参考点也是一种稳健策略
        auto props = faces[fid]->computeProps(nodes);
        x_ref += props.centroid;
    }
    if (!face_indices.empty()) {
        x_ref /= face_indices.size();
    }

    // 2. 累加金字塔
    for (int fid : face_indices) {
        // 获取面的属性
        auto face_props = faces[fid]->computeProps(nodes);

        // 金字塔高度 (投影到法向)
        // 注意：假设面法向是朝外的 (Outward normal)
        double h = (face_props.centroid - x_ref).dot(face_props.normal);

        // 金字塔体积 = 1/3 * BaseArea * Height
        double pyr_vol = (1.0 / 3.0) * face_props.area * h;

        // 金字塔质心 = 3/4 * FaceCentroid + 1/4 * Apex
        Eigen::Vector3d pyr_centroid = 0.75 * face_props.centroid + 0.25 * x_ref;

        // 累加
        volume += pyr_vol;
        moment_1st += pyr_centroid * pyr_vol;
    }

    if (std::abs(volume) > 1e-12) {
        centroid = moment_1st / volume;
    } else {
        centroid = x_ref;
    }
}
void VEMMesh::buildTopology()
{
    // 1. 重置并预分配
    face_to_elements_.clear();
    face_to_elements_.resize(faces_.size());

    // 2. 遍历所有单元，填充 Face -> Element 映射
    for (const auto& elem : elements_) {
        // 如果单元被标记为无效（未来做合并时会用到），则跳过
        // 目前没有 active 标志，暂不处理

        for (int fid : elem.face_indices) {
            if (fid < 0 || fid >= faces_.size()) {
                std::cerr << "[Error] Element " << elem.id << " refers to invalid face " << fid << "\n";
                continue;
            }
            face_to_elements_[fid].push_back(elem.id);
        }
    }
}

const std::vector<int>& VEMMesh::getElementsOnFace(int face_id) const
{
    if (face_id < 0 || face_id >= face_to_elements_.size()) {
        static std::vector<int> empty;
        return empty;
    }
    return face_to_elements_[face_id];
}

std::vector<int> VEMMesh::getElementNeighbors(int elem_id) const
{
    std::vector<int> neighbors;
    if (elem_id < 0 || elem_id >= elements_.size()) {
        return neighbors;
    }

    const auto& elem = elements_[elem_id];

    // 遍历该单元的所有面
    for (int fid : elem.face_indices) {
        const auto& shared_elems = getElementsOnFace(fid);

        // 查找共享该面的其他单元
        for (int neighbor_id : shared_elems) {
            if (neighbor_id != elem_id) {
                neighbors.push_back(neighbor_id);
            }
        }
    }

    // 去重 (因为可能通过多个面连接到同一个邻居? 对于凸多面体一般不会，但为了健壮性)
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

    return neighbors;
}

bool VEMMesh::checkManifold() const
{
    int boundary_count = 0;
    int internal_count = 0;
    int error_count = 0;

    for (size_t fid = 0; fid < face_to_elements_.size(); ++fid) {
        size_t count = face_to_elements_[fid].size();
        if (count == 1) {
            boundary_count++;
        } else if (count == 2) {
            internal_count++;
        } else if (count > 2) {
            error_count++;
            // std::cerr << "[Topology Error] Face " << fid << " is shared by " << count << " elements!" << std::endl;
        } else {
            // count == 0: 孤立的面，未被任何单元使用
        }
    }

    std::cout << "[Topology Check] Boundary Faces: " << boundary_count
              << ", Internal Faces: " << internal_count
              << ", Non-Manifold/Error: " << error_count << "\n";

    return (error_count == 0);
}
int VEMMesh::addPolygonFace(const std::vector<int>& ids)
{
    std::vector<int> sorted_nodes = ids;
    std::sort(sorted_nodes.begin(), sorted_nodes.end());

    auto it = face_map_.find(sorted_nodes);
    if (it != face_map_.end()) {
        return it->second;
    }

    int new_id = faces_.size();

    auto face = std::make_unique<PolygonFace>(ids);
    faces_.push_back(std::move(face));

    face_map_[sorted_nodes] = new_id;

    return new_id;
}

int VEMMesh::addTriFace(const std::vector<int>& ids)
{
    std::vector<int> sorted_nodes = ids;
    std::sort(sorted_nodes.begin(), sorted_nodes.end());

    auto it = face_map_.find(sorted_nodes);
    if (it != face_map_.end()) {
        return it->second;
    }

    faces_.push_back(std::make_unique<TriFace>(ids));

    face_map_[sorted_nodes] = faces_.size() - 1;
    return faces_.size() - 1;
}
