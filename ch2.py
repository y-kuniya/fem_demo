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


class Mesh:
    def __init__(self):
        self.nodes_dict = {}
        self.length = None
        self.height = None
        self.mesh_size = None
        self.num_mesh_len = None
        self.num_mesh_hei = None

    def generate_nodes(self, mesh_size, length, height):
        self.mesh_size = mesh_size
        self.length = length
        self.height = height

        self.num_mesh_len = int(length / mesh_size)
        self.num_mesh_hei = int(height / mesh_size)

        if self.num_mesh_len == 0:
            self.num_mesh_len = 1
        if self.num_mesh_hei == 0:
            self.num_mesh_hei = 1

        x = np.linspace(0, length, self.num_mesh_len + 1)
        y = np.linspace(0, height, self.num_mesh_hei + 1)
        X, Y = np.meshgrid(x, y)
        X = X.ravel()
        Y = Y.ravel()

        self.nodes_dict = {
            i: Node(i, x_coord, y_coord)
            for i, (x_coord, y_coord) in enumerate(zip(X, Y))
        }
        
        for node in self.nodes_dict.values():
            print(node)

mesh = Mesh()
mesh.generate_nodes(mesh_size=5., length=15., height=15.)