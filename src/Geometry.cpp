#include "Geometry.hpp"
#include <cmath>
#include <iostream>
#include <numeric>

// ================= VEMFace Implementations =================

std::unique_ptr<VEMFace> TriFace::clone() const
{
    // 调用拷贝构造函数创建新对象，并封装进 unique_ptr
    return std::make_unique<TriFace>(*this);
}

GeometricProps TriFace::computeProps(const Eigen::MatrixXd& all_nodes) const
{
    Eigen::Vector3d p0 = all_nodes.col(node_indices[0]);
    Eigen::Vector3d p1 = all_nodes.col(node_indices[1]);
    Eigen::Vector3d p2 = all_nodes.col(node_indices[2]);

    Eigen::Vector3d v1 = p1 - p0;
    Eigen::Vector3d v2 = p2 - p0;
    Eigen::Vector3d normal = v1.cross(v2);
    double area = 0.5 * normal.norm();

    if (area > 1e-12)
        normal /= (2.0 * area);
    else
        normal.setZero();

    Eigen::Vector3d centroid = (p0 + p1 + p2) / 3.0;

    return { area, centroid, normal };
}

Eigen::Vector4d TriFace::integrateMonomials(const Eigen::MatrixXd& all_nodes) const
{
    auto props = computeProps(all_nodes);
    double area = props.area;
    Eigen::Vector3d c = props.centroid;

    return Eigen::Vector4d(area, area * c.x(), area * c.y(), area * c.z());
}

std::map<int, double> TriFace::computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const
{
    std::map<int, double> res;
    auto props = computeProps(all_nodes);
    double area = props.area;

    for (int idx : node_indices) {
        res[idx] = area / 3.0;
    }
    return res;
}

// ---------------------------------------------------------

std::unique_ptr<VEMFace> PolygonFace::clone() const
{
    return std::make_unique<PolygonFace>(*this);
}

GeometricProps PolygonFace::computeProps(const Eigen::MatrixXd& all_nodes) const
{
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int idx : node_indices)
        center += all_nodes.col(idx);
    center /= node_indices.size();

    double total_area = 0.0;
    Eigen::Vector3d weighted_normal = Eigen::Vector3d::Zero();
    Eigen::Vector3d weighted_centroid = Eigen::Vector3d::Zero();

    int n = node_indices.size();
    for (int i = 0; i < n; ++i) {
        Eigen::Vector3d p1 = all_nodes.col(node_indices[i]);
        Eigen::Vector3d p2 = all_nodes.col(node_indices[(i + 1) % n]);

        Eigen::Vector3d v1 = p1 - center;
        Eigen::Vector3d v2 = p2 - center;
        Eigen::Vector3d sub_normal = v1.cross(v2);
        double sub_area = 0.5 * sub_normal.norm();

        total_area += sub_area;
        weighted_normal += sub_normal;

        Eigen::Vector3d sub_centroid = (center + p1 + p2) / 3.0;
        weighted_centroid += sub_centroid * sub_area;
    }

    Eigen::Vector3d normal = weighted_normal.normalized();
    Eigen::Vector3d centroid = (total_area > 1e-12) ? (weighted_centroid / total_area) : center;

    return { total_area, centroid, normal };
}

Eigen::Vector4d PolygonFace::integrateMonomials(const Eigen::MatrixXd& all_nodes) const
{
    auto props = computeProps(all_nodes);
    return Eigen::Vector4d(props.area,
        props.area * props.centroid.x(),
        props.area * props.centroid.y(),
        props.area * props.centroid.z());
}

std::map<int, double> PolygonFace::computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const
{
    std::map<int, double> res;
    for (int idx : node_indices)
        res[idx] = 0.0;

    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    for (int idx : node_indices)
        center += all_nodes.col(idx);
    center /= node_indices.size();

    int n = node_indices.size();
    for (int i = 0; i < n; ++i) {
        int idx1 = node_indices[i];
        int idx2 = node_indices[(i + 1) % n];

        Eigen::Vector3d p1 = all_nodes.col(idx1);
        Eigen::Vector3d p2 = all_nodes.col(idx2);

        double sub_area = 0.5 * (p1 - center).cross(p2 - center).norm();
        double contribution_from_center = (sub_area / 3.0) / n;

        res[idx1] += sub_area / 3.0;
        res[idx2] += sub_area / 3.0;

        for (int k : node_indices) {
            res[k] += contribution_from_center;
        }
    }

    return res;
}

// ================= Geometry Utilities Implementation =================

namespace GeometryUtils {

void computePolyhedronProps(
    const Eigen::MatrixXd& all_nodes,
    const std::vector<const VEMFace*>& faces,
    double& out_volume,
    Eigen::Vector3d& out_centroid)
{
    double total_surf_area = 0.0;
    Eigen::Vector3d surf_centroid = Eigen::Vector3d::Zero();

    for (const auto* face : faces) {
        if (!face)
            continue;
        auto props = face->computeProps(all_nodes);
        surf_centroid += props.centroid * props.area;
        total_surf_area += props.area;
    }

    if (total_surf_area > 1e-12) {
        surf_centroid /= total_surf_area;
    } else {
        surf_centroid = Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d x0 = surf_centroid;

    double total_volume = 0.0;
    Eigen::Vector3d moment = Eigen::Vector3d::Zero();

    for (const auto* face : faces) {
        if (!face)
            continue;
        auto props = face->computeProps(all_nodes);

        Eigen::Vector3d outward_normal = props.normal;
        if ((props.centroid - x0).dot(outward_normal) < 0) {
            outward_normal = -outward_normal;
        }

        double h = (props.centroid - x0).dot(outward_normal);
        double pyr_vol = (1.0 / 3.0) * props.area * h;
        Eigen::Vector3d pyr_centroid = 0.75 * props.centroid + 0.25 * x0;

        total_volume += pyr_vol;
        moment += pyr_centroid * pyr_vol;
    }

    out_volume = total_volume;
    if (std::abs(total_volume) > 1e-12) {
        out_centroid = moment / total_volume;
    } else {
        out_centroid = x0;
    }
}

} // namespace GeometryUtils
