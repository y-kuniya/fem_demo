import pygmsh
import meshio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

with pygmsh.geo.Geometry() as geom:
    # 台形の頂点を指定
    points = [
        [0.0, 0.0, 0.0],    # 左下の頂点
        [2.0, 0.0, 0.0],    # 右下の頂点
        [1.5, 1.0, 0.0],    # 右上の頂点
        [0.5, 1.0, 0.0],    # 左上の頂点
    ]
    
    # メッシュサイズを設定（荒めに）
    mesh_size = 0.2
    
    # 台形のポリゴンを作成
    polygon = geom.add_polygon(points, mesh_size=mesh_size)
    
    # メッシュを生成
    mesh = geom.generate_mesh()

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
plt.scatter(points[:, 0], points[:, 1], color='blue', s=40, label='Nodes')

# 各ノードに番号を表示
for i, (x, y) in enumerate(points):
    plt.text(x, y, str(i), color='red', fontsize=10, ha='center', va='center')

# 台形の輪郭を強調表示
trapezoid_outline = np.array([[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.5, 1.0], [0.0, 0.0]])
plt.plot(trapezoid_outline[:, 0], trapezoid_outline[:, 1], 'k-', linewidth=2)

# プロットの設定
plt.gca().set_aspect('equal')
plt.title("Trapezoid Mesh with All Nodes")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig("./data/trapezoid_mesh.png", dpi=300)
plt.show()