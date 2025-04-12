from triangle_mesh import Mesh
import numpy as np 

h = 5.
mesh = Mesh()
mesh.setup_mesh(mesh_size=5., length=10., height=10.)
mesh.get_element_coord()
# mesh.plot_mesh()

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
print("A")
print(A) 
print("f")
print(f)
            
