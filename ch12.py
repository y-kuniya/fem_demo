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

# 2D座標と三角形情報を取得
points = mesh.points[:, :2]
elements = mesh.cells_dict["triangle"]

# 境界線のノードを取得
boundary_nodes = set()
for i in range(1, len(lines) + 1):
    # 各境界線に対応するノードを抽出
    boundary_lines = mesh.cell_sets_dict[f"boundary_{i}"]["line"]
    boundary_nodes.update(np.unique(boundary_lines.flatten()))

# 境界ノードのインデックス
boundary_nodes = np.array(list(boundary_nodes))

print("境界ノード:", boundary_nodes)



import matplotlib.pyplot as plt
# 2D座標と三角形情報を取得
points = mesh.points[:, :2]  # Z座標は無視して2Dにする
triangles = mesh.cells_dict["triangle"]

print(f"ノード数: {len(points)}")
print(f"三角形要素数: {len(triangles)}")
# 描画
plt.figure(figsize=(10, 8))

# すべての三角形を表示
plt.triplot(points[:, 0], points[:, 1], triangles, color='gray', alpha=0.6)

# すべてのノードを表示
plt.scatter(points[:, 0], points[:, 1], color='blue', s=120, label='Nodes')

# 各ノードに番号を表示
for i, (x, y) in enumerate(points):
    plt.text(x, y, str(i), color='white', fontsize=10, ha='center', va='center')

# プロットの設定
plt.gca().set_aspect('equal')
plt.title("Trapezoid Mesh with All Nodes")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("./data/trapezoid_mesh.png", dpi=300)
plt.show()