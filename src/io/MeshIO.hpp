#pragma once
#include <string>
#include <memory>
#include "../VEMMesh.hpp"

namespace vem::io {

// 抽象接口
class MeshIO {
public:
    virtual ~MeshIO() = default;
    
    // 加载网格
    virtual std::unique_ptr<VEMMesh> load(const std::string& filename) const = 0;
    
    // 写入网格
    virtual void save(const VEMMesh& mesh, const std::string& filename) const = 0;
};

// 工厂方法 (可选，暂时可以不用，直接实例化具体类)

} // namespace vem::io