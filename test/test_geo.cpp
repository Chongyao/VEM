#include <iostream>
#include <cassert>
#include "io/VTKWriter.hpp"

bool isApprox(double a, double b, double tol = 1e-10) {
    return std::abs(a - b) < tol;
}

bool isApproxVec(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double tol = 1e-10) {
    return (a - b).norm() < tol;
}

int main() {
    VEMMesh mesh;

    // 1. 创建节点 (单位立方体)
    // 0~3: Bottom z=0, 4~7: Top z=1
    Eigen::MatrixXd V(3, 8);
    V << 0, 1, 1, 0, 0, 1, 1, 0,
         0, 0, 1, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 1, 1, 1;
    mesh.setNodes(V);

    // 2. 创建面 (6个四边形面)，注意顺序要保证法向朝外 (CCW viewing from outside)
    // Bottom (z=0): 0-3-2-1 (Normal 0,0,-1)
    // Top (z=1): 4-5-6-7 (Normal 0,0,1)
    // Front (y=0): 0-1-5-4 (Normal 0,-1,0)
    // Right (x=1): 1-2-6-5 (Normal 1,0,0)
    // Back (y=1): 2-3-7-6 (Normal 0,1,0)
    // Left (x=0): 3-0-4-7 (Normal -1,0,0)

    std::vector<int> f1 = {0, 3, 2, 1}; 
    std::vector<int> f2 = {4, 5, 6, 7}; 
    std::vector<int> f3 = {0, 1, 5, 4}; 
    std::vector<int> f4 = {1, 2, 6, 5}; 
    std::vector<int> f5 = {2, 3, 7, 6}; 
    std::vector<int> f6 = {3, 0, 4, 7}; 

    std::vector<int> face_ids;
    face_ids.push_back(mesh.addPolygonFace(f1));
    face_ids.push_back(mesh.addPolygonFace(f2));
    face_ids.push_back(mesh.addPolygonFace(f3));
    face_ids.push_back(mesh.addPolygonFace(f4));
    face_ids.push_back(mesh.addPolygonFace(f5));
    face_ids.push_back(mesh.addPolygonFace(f6));

    // 3. 创建单元 (这会自动触发 updateGeometricProps)
    mesh.addElement(face_ids);

    // 4. 获取单元并验证
    const auto& elem = mesh.getElement(0);
    
    std::cout << "Computed Volume: " << elem.volume << std::endl;
    std::cout << "Computed Centroid: " << elem.centroid.transpose() << std::endl;

    // 验证
    if (isApprox(elem.volume, 1.0) && 
        isApproxVec(elem.centroid, Eigen::Vector3d(0.5, 0.5, 0.5))) {
        std::cout << "[PASSED] Polyhedron Volume/Centroid test passed!" << std::endl;
    } else {
        std::cerr << "[FAILED] Polyhedron geometry mismatch!" << std::endl;
        std::cerr << "Expected Volume: 1.0, Centroid: 0.5 0.5 0.5" << std::endl;
        return 1;
    }
    
    // 5. 导出 VTK 检查 (可选)
    vem::io::VTKWriter writer;
    writer.save(mesh, "cube_test.vtu");

    return 0;
}