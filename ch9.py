import pygmsh
import meshio
import numpy as np 

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
elements = mesh.cells_dict["triangle"]

dim = len(mesh.points)
A = np.zeros((dim, dim))
f = np.zeros(dim)

# 1. 係数行列、係数ベクトルの計算
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


# 2. ディリクレ境界条件
tol = 1e-10  # 浮動小数点誤差を考慮したトレランス
dirichlet_nodes = []

# for i, (x, y) in enumerate(points):
#     if abs(x - 0.0) < tol or abs(y - 0.0) < tol:
#         dirichlet_nodes.append(i)
for i, (x, y) in enumerate(points):
    if abs(x - 0.0) < tol or abs(x - 1.0) < tol or abs(y - 0.0) < tol or abs(y - 1.0) < tol:
        dirichlet_nodes.append(i)

dirichlet_nodes = np.array(dirichlet_nodes)
print(dirichlet_nodes)

# 3. 自由度ノード
all_nodes = np.arange(dim)                              
free_nodes = np.setdiff1d(all_nodes, dirichlet_nodes) 

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


