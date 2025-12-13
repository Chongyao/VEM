#include "VTKWriter.hpp"
#include <fstream>
#include <iostream>
#include <vector>

namespace vem::io {

void VTKWriter::save(const VEMMesh& mesh, const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "[VTKWriter] Error: Could not open " << filename << std::endl;
        return;
    }

    const auto& nodes = mesh.getNodes();
    int n_nodes = mesh.getNumNodes();
    int n_elems = mesh.getNumElements();

    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_elems << "\">\n";

    // Points
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < n_nodes; ++i) {
        out << nodes(0, i) << " " << nodes(1, i) << " " << nodes(2, i) << " ";
    }
    out << "\n        </DataArray>\n";
    out << "      </Points>\n";

    // Cells
    std::vector<int> connectivity; 
    std::vector<int> offsets;       
    std::vector<int> types;         
    std::vector<int> faces_vec;         
    std::vector<int> faceoffsets;   
    
    int current_offset = 0;
    int current_face_offset = 0;

    for (int i = 0; i < n_elems; ++i) {
        const auto& elem = mesh.getElement(i);
        types.push_back(42); // VTK_POLYHEDRON
        
        int cell_face_data_size = 1; 
        faces_vec.push_back(elem.face_indices.size()); 

        std::vector<int> cell_node_ids;

        for (int fid : elem.face_indices) {
            const auto* face = mesh.getFace(fid);
            const auto& node_ids = face->node_indices;
            
            faces_vec.push_back(node_ids.size());      
            for (int nid : node_ids) {
                faces_vec.push_back(nid);              
                cell_node_ids.push_back(nid);      
            }
            cell_face_data_size += 1 + node_ids.size();
        }
        
        current_face_offset += cell_face_data_size;
        faceoffsets.push_back(current_face_offset);

        for(int nid : cell_node_ids) {
            connectivity.push_back(nid);
        }
        current_offset += cell_node_ids.size();
        offsets.push_back(current_offset);
    }

    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int v : connectivity) out << v << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int v : offsets) out << v << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int v : types) out << v << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"faces\" format=\"ascii\">\n";
    for (int v : faces_vec) out << v << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"Int32\" Name=\"faceoffsets\" format=\"ascii\">\n";
    for (int v : faceoffsets) out << v << " ";
    out << "\n        </DataArray>\n";
    out << "      </Cells>\n";

    // 调试数据：体积
    out << "      <CellData Scalars=\"Volume\">\n";
    out << "        <DataArray type=\"Float64\" Name=\"Volume\" format=\"ascii\">\n";
    for (int i = 0; i < n_elems; ++i) {
        out << mesh.getElement(i).volume << " ";
    }
    out << "\n        </DataArray>\n";
    out << "      </CellData>\n";

    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}

} // namespace vem::io