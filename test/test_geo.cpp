#include <iostream>
#include "VEMMesh.h"

int main() {
    VEMMesh mesh;

    // 1. 创建节点 (单位立方体)
    Eigen::MatrixXd V(3, 8);
    V << 0, 1, 1, 0, 0, 1, 1, 0,
         0, 0, 1, 1, 0, 0, 1, 1,
         0, 0, 0, 0, 1, 1, 1, 1;
    mesh.setNodes(V);

    // 2. 创建面 (6个四边形面)
    // 逆时针顺序定义 (从外面看)
    std::vector<int> f1 = {0, 3, 2, 1}; // Bottom
    std::vector<int> f2 = {4, 5, 6, 7}; // Top
    std::vector<int> f3 = {0, 1, 5, 4}; // Front
    std::vector<int> f4 = {1, 2, 6, 5}; // Right
    std::vector<int> f5 = {2, 3, 7, 6}; // Back
    std::vector<int> f6 = {3, 0, 4, 7}; // Left

    // 这里我们用 addPolygonFace 演示通用性，也可以扩展 addQuadFace
    std::vector<int> face_ids;
    face_ids.push_back(mesh.addPolygonFace(f1));
    face_ids.push_back(mesh.addPolygonFace(f2));
    face_ids.push_back(mesh.addPolygonFace(f3));
    face_ids.push_back(mesh.addPolygonFace(f4));
    face_ids.push_back(mesh.addPolygonFace(f5));
    face_ids.push_back(mesh.addPolygonFace(f6));

    // 3. 创建单元 (1个六面体)
    mesh.addElement(face_ids);

    // 4. 测试几何计算
    auto props = mesh.getFace(0)->computeProps(mesh.getNodes());
    std::cout << "Bottom Face Area: " << props.area << " (Expected: 1.0)" << std::endl;
    std::cout << "Bottom Face Normal: " << props.normal.transpose() << " (Expected: 0 0 -1)" << std::endl;

    // 5. 导出 VTK
    mesh.exportToVTK("cube.vtu");

    return 0;
}