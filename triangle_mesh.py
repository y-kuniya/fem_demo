import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

class Node:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y

    def __repr__(self):
        return f"Node(id={self.id}, x={self.x}, y={self.y})"


class Element:
    def __init__(self, id, node_ids):
        self.id = id
        self.node_ids = node_ids
        self.xy = None
        self.triangle_type = None

    def __repr__(self):
        return f"Element(id={self.id}, node_ids={self.node_ids})"

    def get_coordination(self, nodes_dict):
        res = []
        for node_idx in self.node_ids:
            res.append([nodes_dict[node_idx].x, nodes_dict[node_idx].y])
        self.xy = np.array(res)
    
    def determine_triangle_type(self):
        if self.xy is None:
            return

        # y座標の平均より大きい点が1つ → 上向き三角形
        # else                      → 下向き三角形
        y_vals = self.xy[:, 1]
        y_mean = np.mean(y_vals)
        above = (y_vals > y_mean).sum()

        if above == 1:
            self.triangle_type = "upper"
        else :
            self.triangle_type = "lower"

class Mesh:
    def __init__(self):
        self.nodes_dict = {}
        self.elements = []

    def generate_nodes(self, mesh_size, length, height):
        num_mesh_len = max(int(length / mesh_size), 1)
        num_mesh_hei = max(int(height / mesh_size), 1)

        x = np.linspace(0, length, num_mesh_len + 1)
        y = np.linspace(0, height, num_mesh_hei + 1)
        X, Y = np.meshgrid(x, y)
        X = X.ravel()
        Y = Y.ravel()

        self.nodes_dict = {
            i: Node(i, x_coord, y_coord)
            for i, (x_coord, y_coord) in enumerate(zip(X, Y))
        }

        # 保存しておくなら必要に応じて下記のように定義もできる
        self._num_mesh_len = num_mesh_len
        self._num_mesh_hei = num_mesh_hei

    def generate_triangular_elements(self):
        num_mesh_len = self._num_mesh_len
        num_mesh_hei = self._num_mesh_hei

        elem_idx = 0
        triangle_elems = []

        for i in range(num_mesh_hei):
            for j in range(num_mesh_len):
                n1 = i * (num_mesh_len + 1) + j
                n2 = n1 + 1
                n3 = n2 + (num_mesh_len + 1)
                n4 = n1 + (num_mesh_len + 1)

                triangle_elems.append(Element(elem_idx, [n1, n2, n3]))
                elem_idx += 1
                triangle_elems.append(Element(elem_idx, [n1, n3, n4]))
                elem_idx += 1

        self.elements = triangle_elems
    
    def setup_mesh(self, mesh_size, length, height):
        self.generate_nodes(mesh_size, length, height)
        self.generate_triangular_elements()
        # for node in self.nodes_dict.values():
        #     print(node)
        # for elem in self.elements:
        #     print(elem)
    
    def get_element_coord(self):
        for elm in self.elements:
            elm.get_coordination(self.nodes_dict)
            elm.determine_triangle_type()
    
    def print_element_coordinates(self):
        for elem in self.elements:
            print(f"Element {elem.id}:")
            print(elem.triangle_type)
            print(elem.xy)

    def plot_mesh(self):
        # ノード座標を numpy 配列に変換
        coords = np.array([[node.x, node.y] for node in self.nodes_dict.values()])
        x = coords[:, 0]
        y = coords[:, 1]

        # 三角形のノードインデックスリストを作成
        triangles = [elem.node_ids for elem in self.elements]

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

        # ノード番号を描画
        for i, (xi, yi) in enumerate(zip(x, y)):
            plt.text(xi, yi, str(i), color='blue', fontsize=12)

        plt.savefig("./data/mesh.png")
