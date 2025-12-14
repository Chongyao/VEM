#pragma once
#include <Eigen/Dense>
#include <map>
#include <memory>
#include <vector>

// Forward declarations
class VEMFace;

// [Moved from VEMMesh.hpp]
// 移动到这里是为了解决循环依赖，并让 GeometryUtils 可以直接操作它
struct PolyhedronElement {
    int id;
    std::vector<int> face_indices;
    bool active = true; // [New] Active flag for CutVEM

    // Cached geometric data
    double volume = 0.0;
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
};

struct GeometricProps {
    double area;
    Eigen::Vector3d centroid;
    Eigen::Vector3d normal;
};

class VEMFace {
public:
    std::vector<int> node_indices;

    VEMFace(const std::vector<int>& nodes)
        : node_indices(nodes)
    {
    }
    virtual ~VEMFace() = default;

    // [Crucial] 纯虚函数 clone，用于支持 VEMMesh 的深拷贝
    virtual std::unique_ptr<VEMFace> clone() const = 0;

    virtual GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const = 0;
    virtual Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const = 0;
    virtual std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const = 0;
};

class TriFace : public VEMFace {
public:
    TriFace(const std::vector<int>& nodes)
        : VEMFace(nodes)
    {
    }

    std::unique_ptr<VEMFace> clone() const override;

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};

class PolygonFace : public VEMFace {
public:
    PolygonFace(const std::vector<int>& nodes)
        : VEMFace(nodes)
    {
    }

    std::unique_ptr<VEMFace> clone() const override;

    GeometricProps computeProps(const Eigen::MatrixXd& all_nodes) const override;
    Eigen::Vector4d integrateMonomials(const Eigen::MatrixXd& all_nodes) const override;
    std::map<int, double> computeShapeFuncIntegrals(const Eigen::MatrixXd& all_nodes) const override;
};

// [New] Geometry Utilities Namespace
// 提供无状态的几何计算函数，供 VEMMesh 和 AgglomerationManager 复用
namespace GeometryUtils {
/**
 * @brief 计算多面体的体积和质心
 * @param all_nodes 全局节点列表
 * @param faces 构成该多面体的面指针列表
 * @param out_volume 输出体积 (引用返回)
 * @param out_centroid 输出质心 (引用返回)
 */
void computePolyhedronProps(
    const Eigen::MatrixXd& all_nodes,
    const std::vector<const VEMFace*>& faces,
    double& out_volume,
    Eigen::Vector3d& out_centroid);
}
