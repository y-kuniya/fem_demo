import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.tri as tri

class Node:
    def __init__(self,idx,x,y):
        self.id = idx 
        self.x  = x
        self.y  = y 
    def __str__(self):
        return 'Node {}: x {}, y {}'.format(self.id, self.x, self.y)
    
class Element:
    def __init__(self, idx, node_idxs):
        self.id = idx
        self.node = node_idxs
    
    def get_coordination(self, nodes_dict):
        res = []
        for node_idx in self.node:
            res.append([nodes_dict[node_idx].x, nodes_dict[node_idx].y])
        self.xy = np.array(res)
    
    def __str__(self):
        return 'Element {}: Nodes {}'.format(self.id, self.node)
    
class Mesh:
    def __init__(self, nodes_dict, elements):
        self.nodes = nodes_dict
        self.elements = elements
        self.get_element_coord()
    
    def get_element_coord(self):
        for elm in self.elements:
            elm.get_coordination(self.nodes)


mesh_size = 5.
length = 20.
height = 20.
# Nodeの作成
num_mesh_len = int(length / mesh_size)
num_mesh_hei = int(height / mesh_size)
if num_mesh_len == 0:
    num_mesh_len = 1
if num_mesh_hei == 0:
    num_mesh_hei = 1

x = np.linspace(0, length, num_mesh_len + 1)
y = np.linspace(0, height, num_mesh_hei + 1)
X, Y = np.meshgrid(x, y)
X = X.ravel()
Y = Y.ravel()

nodes_dict = {}
for i, coord in enumerate(zip(X, Y)):
    nodes_dict[i] = Node(i, coord[0], coord[1])

nodes_dict = {}
for i, coord in enumerate(zip(X, Y)):
    nodes_dict[i] = Node(i, coord[0], coord[1])

for node in nodes_dict.values():
    print(node)

triangle_elems = []
elem_idx = 0

for i in range(num_mesh_hei):
    for j in range(num_mesh_len):
        n1 = i * (num_mesh_len + 1) + j
        n2 = n1 + 1
        n3 = n2 + (num_mesh_len + 1)
        n4 = n1 + (num_mesh_len + 1)

        # 三角形1（上半分）
        triangle_elems.append(Element(elem_idx, [n1, n2, n3]))
        elem_idx += 1

        # 三角形2（下半分）
        triangle_elems.append(Element(elem_idx, [n1, n3, n4]))
        elem_idx += 1

for elem in triangle_elems:
    print(elem)

# ノード座標を numpy 配列に変換
coords = np.array([[node.x, node.y] for node in nodes_dict.values()])
x = coords[:, 0]
y = coords[:, 1]

# 三角形のノードインデックスリストを作成
triangles = [elem.node for elem in triangle_elems]

# Triangulation を作成
triang = tri.Triangulation(x, y, triangles)

# プロット
plt.figure(figsize=(6, 6))
plt.triplot(triang, color='black')
plt.gca().set_aspect('equal')
plt.title("Triangular Mesh")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
for i, (xi, yi) in enumerate(zip(x, y)):
    plt.text(xi, yi, str(i), color='blue', fontsize=12)
plt.savefig("mesh.png")