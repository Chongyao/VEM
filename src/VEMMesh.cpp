#include "VEMMesh.hpp"
#include <algorithm>
#include <iostream>
#include <set>

// [New] Copy Constructor - Deep Copy
VEMMesh::VEMMesh(const VEMMesh& other)
    : nodes_(other.nodes_)
    , elements_(other.elements_)
    , // vector copy (copies active flags etc.)
    face_to_elements_(other.face_to_elements_)
    , face_map_(other.face_map_)
{
    // unique_ptr 不能拷贝，必须手动克隆
    faces_.reserve(other.faces_.size());
    for (const auto& face : other.faces_) {
        faces_.push_back(face->clone());
    }
}

// [New] Copy Assignment - Deep Copy
VEMMesh& VEMMesh::operator=(const VEMMesh& other)
{
    if (this != &other) {
        nodes_ = other.nodes_;
        elements_ = other.elements_;
        face_to_elements_ = other.face_to_elements_;
        face_map_ = other.face_map_;

        // Clear old faces and clone new ones
        faces_.clear();
        faces_.reserve(other.faces_.size());
        for (const auto& face : other.faces_) {
            faces_.push_back(face->clone());
        }
    }
    return *this;
}

int VEMMesh::addTriFace(const std::vector<int>& ids)
{
    return addPolygonFace(ids);
}

int VEMMesh::addPolygonFace(const std::vector<int>& node_indices)
{
    std::vector<int> sorted_nodes = node_indices;
    std::sort(sorted_nodes.begin(), sorted_nodes.end());

    auto it = face_map_.find(sorted_nodes);
    if (it != face_map_.end()) {
        return it->second;
    }

    int new_id = faces_.size();
    if (node_indices.size() == 3) {
        faces_.push_back(std::make_unique<TriFace>(node_indices));
    } else {
        faces_.push_back(std::make_unique<PolygonFace>(node_indices));
    }
    // 确保 node_indices 设置正确
    faces_.back()->node_indices = node_indices;

    face_map_[sorted_nodes] = new_id;
    return new_id;
}

int VEMMesh::addElement(const std::vector<int>& face_ids)
{
    PolyhedronElement elem;
    elem.id = elements_.size();
    elem.face_indices = face_ids;
    elem.active = true;

    // 1. 准备数据调用几何计算
    std::vector<const VEMFace*> face_ptrs;
    face_ptrs.reserve(face_ids.size());
    for (int fid : face_ids) {
        if (fid >= 0 && fid < faces_.size()) {
            face_ptrs.push_back(faces_[fid].get());
        }
    }

    // 2. 自动计算体积和质心
    GeometryUtils::computePolyhedronProps(nodes_, face_ptrs, elem.volume, elem.centroid);

    // 3. 存入列表
    elements_.push_back(elem);
    int new_elem_id = elem.id;

    // 4. 自动更新拓扑
    if (face_to_elements_.size() < faces_.size()) {
        face_to_elements_.resize(faces_.size());
    }
    for (int fid : face_ids) {
        if (fid < face_to_elements_.size()) {
            face_to_elements_[fid].push_back(new_elem_id);
        }
    }

    return new_elem_id;
}

void VEMMesh::deactivateElement(int elem_id)
{
    if (elem_id >= 0 && elem_id < elements_.size()) {
        elements_[elem_id].active = false;
    }
}

bool VEMMesh::isElementActive(int elem_id) const
{
    if (elem_id >= 0 && elem_id < elements_.size()) {
        return elements_[elem_id].active;
    }
    return false;
}

void VEMMesh::buildTopology()
{
    face_to_elements_.clear();
    face_to_elements_.resize(faces_.size());

    for (const auto& elem : elements_) {
        if (!elem.active)
            continue;

        for (int fid : elem.face_indices) {
            if (fid < 0 || fid >= faces_.size()) {
                std::cerr << "[Error] Element " << elem.id << " refers to invalid face " << fid << std::endl;
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
    if (elem_id < 0 || elem_id >= elements_.size())
        return neighbors;
    if (!elements_[elem_id].active)
        return neighbors;

    const auto& elem = elements_[elem_id];

    for (int fid : elem.face_indices) {
        const auto& shared_elems = getElementsOnFace(fid);
        for (int neighbor_id : shared_elems) {
            if (neighbor_id != elem_id && isElementActive(neighbor_id)) {
                neighbors.push_back(neighbor_id);
            }
        }
    }

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
        int active_users = 0;
        for (int eid : face_to_elements_[fid]) {
            if (isElementActive(eid))
                active_users++;
        }

        if (active_users == 1)
            boundary_count++;
        else if (active_users == 2)
            internal_count++;
        else if (active_users > 2)
            error_count++;
    }

    std::cout << "[Topology Check] Boundary Faces: " << boundary_count
              << ", Internal Faces: " << internal_count
              << ", Non-Manifold/Error: " << error_count << std::endl;

    return (error_count == 0);
}
