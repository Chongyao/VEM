import numpy as np
from scipy.spatial import Voronoi
from tqdm import tqdm
from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkCommonDataModel import VTK_POLYHEDRON, vtkUnstructuredGrid
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridWriter
import sys

# ================= 参数 =================
OUTPUT_FILE = "beam_cvt_final.vtu"
N_SEEDS = 2000 # 为了测试稳定，先用少一点，或者用你的 1000
if len(sys.argv) > 1:
    N_SEEDS = int(sys.argv[1])
N_ITERS = 10
BOUNDS = [-2.0, 2.0, -0.5, 0.5, -0.5, 0.5]

def reflect_points(points, bounds):
    """镜像生成幽灵点"""
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    planes = [(0, x_min), (0, x_max), (1, y_min), (1, y_max), (2, z_min), (2, z_max)]
    all_points = [points]
    for axis, val in planes:
        p_copy = points.copy()
        p_copy[:, axis] = 2 * val - p_copy[:, axis]
        all_points.append(p_copy)
    return np.vstack(all_points)

print("--- Generating Normal-Corrected Polyhedral Mesh ---")

# 1. 生成种子
bx_min, bx_max, by_min, by_max, bz_min, bz_max = BOUNDS
seeds = np.random.rand(N_SEEDS, 3)
seeds[:, 0] = seeds[:, 0] * (bx_max - bx_min) + bx_min
seeds[:, 1] = seeds[:, 1] * (by_max - by_min) + by_min
seeds[:, 2] = seeds[:, 2] * (bz_max - bz_min) + bz_min

# Lloyd 迭代
for it in range(N_ITERS):
    pts_mirror = reflect_points(seeds, BOUNDS)
    vor = Voronoi(pts_mirror)
    new_seeds = []
    for i in range(N_SEEDS):
        region = vor.regions[vor.point_region[i]]
        if not region: 
            new_seeds.append(seeds[i])
            continue
        verts = vor.vertices[region]
        center = np.mean(verts, axis=0)
        center[0] = np.clip(center[0], bx_min, bx_max)
        center[1] = np.clip(center[1], by_min, by_max)
        center[2] = np.clip(center[2], bz_min, bz_max)
        new_seeds.append(center)
    seeds = np.array(new_seeds)

# 2. 构建拓扑
print("Constructing Topology...")
pts_mirror = reflect_points(seeds, BOUNDS)
vor = Voronoi(pts_mirror)

cells_faces = {} # seed_idx -> list of faces (each face is a list of node ids)

# 遍历所有脊
for i in range(len(vor.ridge_points)):
    seed_indices = vor.ridge_points[i]      # [id1, id2]
    face_vert_indices = vor.ridge_vertices[i] # [v1, v2, ...]
    
    if -1 in face_vert_indices or len(face_vert_indices) < 3:
        continue

    face_pts = vor.vertices[face_vert_indices]
    
    # --- 排序逻辑 ---
    # 计算脊的基准法向 (从 S1 指向 S2)
    p1 = pts_mirror[seed_indices[0]]
    p2 = pts_mirror[seed_indices[1]]
    normal_ref = p2 - p1
    normal_ref /= np.linalg.norm(normal_ref) + 1e-12
    
    # 对顶点进行逆时针排序 (相对于 normal_ref)
    center = np.mean(face_pts, axis=0)
    axis_u = face_pts[0] - center
    axis_u /= np.linalg.norm(axis_u) + 1e-12
    axis_v = np.cross(normal_ref, axis_u)
    axis_v /= np.linalg.norm(axis_v) + 1e-12
    
    vecs = face_pts - center
    u_proj = np.dot(vecs, axis_u)
    v_proj = np.dot(vecs, axis_v)
    angles = np.arctan2(v_proj, u_proj)
    sorted_order = np.argsort(angles)
    
    # 得到标准顺序 (对于 S1 来说是 Outward 的)
    sorted_indices = [face_vert_indices[k] for k in sorted_order]
    
    s1, s2 = seed_indices
    
    # 分配给 Seed 1 (S1)
    # normal_ref = S2 - S1. S1在中心，S2在外部。
    # 实际上 Voronoi 定义 ridge 垂直于 S1-S2。
    # 从 S1 看去，S2 方向是向外的。
    # 所以 sorted_indices (CCW) 对应的法向是 normal_ref (Outward)。正确。
    if s1 < N_SEEDS:
        if s1 not in cells_faces: cells_faces[s1] = []
        cells_faces[s1].append(sorted_indices)
        
    # 分配给 Seed 2 (S2)
    # 对于 S2，Outward 方向是 S1 - S2 = -normal_ref。
    # 所以我们需要一个顺序，使得其法向是 -normal_ref。
    # 这等价于将 sorted_indices 逆序。
    if s2 < N_SEEDS:
        if s2 not in cells_faces: cells_faces[s2] = []
        cells_faces[s2].append(list(reversed(sorted_indices)))

# 3. 写入 VTK
ugrid = vtkUnstructuredGrid()
vtk_pts = vtkPoints()
node_map = {} 
current_pt_id = 0

for seed_id in tqdm(range(N_SEEDS)):
    if seed_id not in cells_faces: continue
    faces = cells_faces[seed_id]
    if not faces: continue
    
    faceIdList = vtkIdList()
    faceIdList.InsertNextId(len(faces)) # nFaces
    
    for face_verts in faces:
        faceIdList.InsertNextId(len(face_verts)) # nPts
        for vid in face_verts:
            if vid not in node_map:
                p = vor.vertices[vid]
                vtk_pts.InsertNextPoint(p[0], p[1], p[2])
                node_map[vid] = current_pt_id
                current_pt_id += 1
            faceIdList.InsertNextId(node_map[vid])
            
    ugrid.InsertNextCell(VTK_POLYHEDRON, faceIdList)

ugrid.SetPoints(vtk_pts)

print(f"Writing {OUTPUT_FILE} (Cells: {ugrid.GetNumberOfCells()})...")
writer = vtkXMLUnstructuredGridWriter()
writer.SetInputData(ugrid)
writer.SetFileName(OUTPUT_FILE)
writer.SetDataModeToAscii()
writer.Update()
print("Done.")
