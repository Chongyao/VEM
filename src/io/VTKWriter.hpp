#pragma once
#include "MeshIO.hpp"
#include <Eigen/Core> // 需要包含 Eigen

namespace vem::io {

class VTKWriter {
public:
    // 仅保存网格
    void save(const VEMMesh& mesh, const std::string& filename) const;

    // [New] 保存网格 + 位移场
    void save(const VEMMesh& mesh, const Eigen::VectorXd& u, const std::string& filename) const;
};

} // namespace vem::io
