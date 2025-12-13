#pragma once
#include "MeshIO.hpp"

namespace vem::io {

class VTKLoader : public MeshIO {
public:
    std::unique_ptr<VEMMesh> load(const std::string& filename) const override;
    
    // Loader 不需要 save
    void save(const VEMMesh& mesh, const std::string& filename) const override {}
};

} // namespace vem::io