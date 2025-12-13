#pragma once
#include "MeshIO.hpp"

namespace vem::io {

class VTKWriter : public MeshIO {
public:
    // Writer 不需要 load
    std::unique_ptr<VEMMesh> load(const std::string& filename) const override { return nullptr; }
    
    void save(const VEMMesh& mesh, const std::string& filename) const override;
};

} // namespace vem::io