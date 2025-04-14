import pygmsh
import numpy as np
import meshio

# 星形の頂点を計算する関数
def star_points(radius_outer, radius_inner, num_points):
    angle_step = np.pi / num_points  # 角度のステップ
    points = []

    for i in range(num_points):
        # 外側の頂点
        angle = i * 2 * np.pi / num_points
        x_outer = radius_outer * np.cos(angle)
        y_outer = radius_outer * np.sin(angle)
        points.append([x_outer, y_outer])

        # 内側の頂点（外側と内側の間）
        angle = (i + 0.5) * 2 * np.pi / num_points
        x_inner = radius_inner * np.cos(angle)
        y_inner = radius_inner * np.sin(angle)
        points.append([x_inner, y_inner])

    return np.array(points)

with pygmsh.geo.Geometry() as geom:
    mesh_size = 0.2

    # 星形の頂点を定義
    radius_outer = 1.0  # 外側の半径
    radius_inner = 0.4  # 内側の半径
    num_points = 5  # 5角形の星（奇数の頂点数が必要）

    points = star_points(radius_outer, radius_inner, num_points)

    # 頂点を "Point" として追加
    pygmsh_points = [geom.add_point(p, mesh_size) for p in points]

    # 頂点を結ぶ線を定義
    lines = []
    for i in range(len(pygmsh_points)):
        lines.append(geom.add_line(pygmsh_points[i], pygmsh_points[(i + 1) % len(pygmsh_points)]))

    # 曲線ループと面の定義
    loop = geom.add_curve_loop(lines)
    surface = geom.add_plane_surface(loop)

    # 物理領域を追加（後で参照する用）
    geom.add_physical(surface, "domain")
    for i, line in enumerate(lines):
        geom.add_physical(line, f"boundary_{i+1}")  # 境界線を物理的に指定

    mesh = geom.generate_mesh()

# 0. 2D座標と三角形情報を取得
points = mesh.points[:, :2]
elements = mesh.cells_dict["triangle"]

# 1. 次元と全体の行列を取得
dim = len(mesh.points)
A = np.zeros((dim, dim))
f = np.zeros(dim) 

# 2. 境界線のノードを取得 (ディリクレ境界条件のため)
boundary_nodes = set()
for i in range(1, len(lines) + 1):
    # 各境界線に対応するノードを抽出
    boundary_lines = mesh.cell_sets_dict[f"boundary_{i}"]["line"]
    boundary_nodes.update(np.unique(boundary_lines.flatten()))

# 境界ノードのインデックス
boundary_nodes = np.array(list(boundary_nodes))

# 2. 自由度ノード
all_nodes = np.arange(dim)                              
free_nodes = np.setdiff1d(all_nodes, boundary_nodes)

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
plt.show()
plt.savefig("heat_map_star.png")
