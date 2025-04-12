from triangle_mesh import Mesh
import numpy as np 

h = 0.1
mesh = Mesh()
mesh.setup_mesh(mesh_size=h, length=1., height=1.)
mesh.get_element_coord()
mesh.plot_mesh()

Aupper = np.array([
        [1, -1, 0],
        [-1, 2, -1],
        [0, -1, 1]
        ],dtype=float)
Aupper *= 0.5 

Alower = np.array([
        [1, 0, -1],
        [0, 1, -1],
        [-1, -1,2]
        ],dtype=float)
Alower *= 0.5 

dim = len(mesh.nodes_dict)
A = np.zeros((dim, dim))
f = np.zeros(dim)

for elem in mesh.elements:
    nodes = elem.node_ids
    type  = elem.triangle_type
    Ae = None
    if (type == "upper"):
        Ae = Aupper
    else:
        Ae = Alower

    for k in range(3):
        for l in range(3):
            A[nodes[k],nodes[l]] += Ae[k,l]

    for k in range(3):
        f[nodes[k]] += h*h/6

# print("A")
# print(A) 
# print("f")
# print(f)

# 1. ディリクレ境界ノード
dirichlet_nodes = mesh.get_boundary_nodes_2()             # => [0 3 8]

# 2. 自由度ノード
all_nodes = np.arange(dim)                              # => [0 1 2 3 4 5 6 7 8]
free_nodes = np.setdiff1d(all_nodes, dirichlet_nodes)   # => [1 2 4 5 6 7]

# 3. 縮約
A_reduced = A[np.ix_(free_nodes, free_nodes)]
f_reduced = f[free_nodes]
# print("A_reduced")
# print(A_reduced)
# print("f_reduced")
# print(f_reduced)

# 4. 解く
u_reduced = np.linalg.solve(A_reduced, f_reduced)

# 5. 埋め戻し
u = np.zeros(dim)
u[free_nodes] = u_reduced
# print("u")
# print(u)
print(u[dim//2])
      
