import pygmsh
import meshio

with pygmsh.geo.Geometry() as geom:
    # 四角形の頂点を指定
    points = [
        [0.0, 0.0, 0.0],  # (xmin, ymin)
        [1.0, 0.0, 0.0],  # (xmax, ymin)
        [1.0, 1.0, 0.0],  # (xmax, ymax)
        [0.0, 1.0, 0.0],  # (xmin, ymax)
    ]
    # ポリゴンを作成
    geom.add_polygon(points, mesh_size=0.1)
    mesh = geom.generate_mesh()

# numpyで処理したいとき
points = mesh.points[:, :2]
cells = mesh.cells_dict["triangle"]

# meshioで保存も可
# mesh.write("mesh.vtk")


import matplotlib.pyplot as plt
import matplotlib.tri as tri

# 2D座標と三角形情報を取得
points = mesh.points[:, :2]
triangles = mesh.cells_dict["triangle"]

# Triangulation を作成
triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)

# 描画して保存
plt.figure(figsize=(6, 6))
plt.triplot(triang, color='gray')
plt.gca().set_aspect('equal')
plt.title("Generated Mesh")
plt.xlabel("x")
plt.ylabel("y")
# グリッド線をグレーに設定
plt.grid(True)
# 各ノードの位置に番号を表示
for i, (x, y) in enumerate(points):
    plt.text(x, y, str(i), color='blue', fontsize=12, ha='center', va='center')
plt.savefig("./data/mesh_from_library.png")

