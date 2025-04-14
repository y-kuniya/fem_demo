import pygmsh
import meshio
import numpy as np

with pygmsh.geo.Geometry() as geom:
    mesh_size = 0.2

    # 台形の頂点を Point として定義
    p0 = geom.add_point([0.0, 0.0, 0.0], mesh_size)
    p1 = geom.add_point([2.0, 0.0, 0.0], mesh_size)
    p2 = geom.add_point([1.5, 1.0, 0.0], mesh_size)
    p3 = geom.add_point([0.5, 1.0, 0.0], mesh_size)

    # 各辺を Line として定義
    l0 = geom.add_line(p0, p1)  # bottom
    l1 = geom.add_line(p1, p2)  # right slant
    l2 = geom.add_line(p2, p3)  # top
    l3 = geom.add_line(p3, p0)  # left slant

    # Curve loop & surface
    loop = geom.add_curve_loop([l0, l1, l2, l3])
    surface = geom.add_plane_surface(loop)

    # 物理領域を追加（後で参照する用）
    geom.add_physical(surface, "domain")
    geom.add_physical(l0, "bottom")
    geom.add_physical(l1, "right")
    geom.add_physical(l2, "top")
    geom.add_physical(l3, "left")

    mesh = geom.generate_mesh()


# 2D座標と三角形情報を取得
points = mesh.points[:, :2]
elements = mesh.cells_dict["triangle"]

top_nodes = np.unique(mesh.cell_sets_dict["top"]["line"].flatten())
bottom_nodes = np.unique(mesh.cell_sets_dict["bottom"]["line"].flatten())
left_nodes = np.unique(mesh.cell_sets_dict["left"]["line"].flatten())
right_nodes = np.unique(mesh.cell_sets_dict["right"]["line"].flatten())

#0. 次元と全体の行列を取得
dim = len(mesh.points)
A = np.zeros((dim, dim))
f = np.zeros(dim)

# 1. 全ディリクレ対象ノードのインデックス
dirichlet_nodes = np.unique(np.concatenate([top_nodes, bottom_nodes, left_nodes, right_nodes]))
print(dirichlet_nodes)

# 2. 自由度ノード
all_nodes = np.arange(dim)                              
free_nodes = np.setdiff1d(all_nodes, dirichlet_nodes) 

# 3. 係数行列、係数ベクトルの計算
for elem in elements:
    x0 = points[elem[0]]
    x1 = points[elem[1]]
    x2 = points[elem[2]]

    d1 = x1 - x0
    d2 = x2 - x0
    area = 0.5 * abs(np.cross(d1, d2))

    # 行列要素の計算
    # 1.  
    P = np.ones((3, 3))
    P[0, 1:] = x0  # x, y of point 0
    P[1, 1:] = x1  # x, y of point 1
    P[2, 1:] = x2  # x, y of point 2

    # 2. 逆行列の計算
    if np.abs(np.linalg.det(P)) < 1e-12:
        raise ValueError("行列Pはほぼ特異行列です。")
    else:
        P_inv = np.linalg.inv(P)
    
    #3. 係数行列の計算 
    v1 = P_inv[1]
    v2 = P_inv[2]
    Ae = (np.outer(v1, v1) + np.outer(v2, v2))*area

    #4. 係数ベクトルの計算
    fe = (area / 3.0) * np.ones(3)

    #5. 最後に全体の配列に格納
    for k in range(3):
        for l in range(3):
            A[elem[k],elem[l]] += Ae[k,l]

    for k in range(3):
        f[elem[k]] += fe[k]

# 4. 縮約
A_reduced = A[np.ix_(free_nodes, free_nodes)]
f_reduced = f[free_nodes]

# 5. 解く
u_reduced = np.linalg.solve(A_reduced, f_reduced)

# 6. 埋め戻し
u = np.zeros(dim)
u[free_nodes] = u_reduced


u_max = np.max(u)
print("uの最大値:", u_max)
i_max = np.argmax(u)
print("最大値の位置 (ノード番号):", i_max)
print("座標:", points[i_max])


# 最後に可視化
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# 1. メッシュ構造を作成
triangles = elements  
triang = mtri.Triangulation(points[:, 0], points[:, 1], triangles)

# 2. 各三角形に対応する u の平均値を計算
u_avg = np.mean(u[triangles], axis=1)  # 三角形ごとに頂点のuを平均

# 3. 描画
fig, ax = plt.subplots(figsize=(6, 5))
tpc = ax.tripcolor(triang, facecolors=u_avg, edgecolors='k', cmap='plasma')
ax.set_aspect('equal')
plt.colorbar(tpc, ax=ax, label="u")
ax.set_title("heat map of u")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("heat_map_trapezoid.png")
