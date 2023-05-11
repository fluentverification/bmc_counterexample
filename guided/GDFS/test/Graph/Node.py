class Node:
    def __init__(self):
        self.var_values = tuple()
        self.initial_state = False
        self.out_edges = {}
        self.index = -1
        self.reachability_probability = 0
        self.in_edges = {}

    def add_out_edge(self, edge):
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.out_edges:
            self.out_edges[edge_tuple] = edge

    def add_in_edge(self, edge):
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.in_edges:
            self.in_edges[edge_tuple] = edge

    def __eq__(self, node) -> bool:
        return self.var_values == node.var_values

    def __lt__(self, node):
        return True