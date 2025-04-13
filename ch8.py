from triangle_mesh import Mesh
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

u_list = []
for m in range(2,11):
    h = 1.0/m
    mesh = Mesh()
    mesh.setup_mesh(mesh_size=h, length=1., height=1.)
    mesh.get_element_coord()

    dim = len(mesh.nodes_dict)
    A = np.zeros((dim, dim))
    f = np.zeros(dim)

    # 1. 係数行列、係数ベクトルの計算
    for elem in mesh.elements:
        # 三角形の面積の計算
        a,b,c = np.array(elem.xy)
        d1 = b - a
        d2 = c - a
        area = 0.5 * abs(np.cross(d1, d2))

        # 行列要素の計算
        # 1.  
        P = np.ones((3, 3))
        P[:, 1:] = np.array(elem.xy)

        # 2. 逆行列の計算
        if np.linalg.det(P) != 0:
            P_inv = np.linalg.inv(P)
        else:
            raise ValueError("行列Pは正則ではないので逆行列を持ちません。")
        
        #3. 係数行列の計算 
        v1 = P_inv[1]
        v2 = P_inv[2]
        Ae = (np.outer(v1, v1) + np.outer(v2, v2))*area

        #4. 係数ベクトルの計算
        fe = (area / 3.0) * np.ones(3)

        #5. 最後に全体の配列に格納
        nodes = elem.node_ids
        for k in range(3):
            for l in range(3):
                A[nodes[k],nodes[l]] += Ae[k,l]

        for k in range(3):
            f[nodes[k]] += fe[k]

    # 2. ディリクレ境界条件
    dirichlet_nodes = mesh.get_boundary_nodes_2()

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
    u_list.append(max(u))


df = pd.DataFrame({'u(0.5,0.5)': u_list},index=np.arange(2, 11))

fig, ax = plt.subplots(figsize=(5, 3))  # サイズ指定（任意）
ax.axis('tight')
ax.axis('off')
#table = ax.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center')
#テーブル生成（行ラベル＝mの値を表示）
table = ax.table(
    cellText=df.values,
    rowLabels=df.index,         # ← ここで m=2~10 を行ラベルに
    colLabels=df.columns,
    loc='center',
    cellLoc='center'
)

# PNGファイルとして保存
plt.savefig("u_center_for_various_m.png", bbox_inches='tight', dpi=300)

