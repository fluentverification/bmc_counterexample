class Graph:
    def __init__(self):
        self.nodes = {}
        self.edges = {}

    def get_node(self, node_tuple):
        if node_tuple in self.nodes:
            return (True, self.nodes[node_tuple])
        return (False, None)
    
    def get_edge(self, edge_tuple):
        if edge_tuple in self.edges:
            return (True, self.edges[edge_tuple])
        return (False, None)
    
    def add_node(self, node):
        if node.var_values not in self.nodes:
            self.nodes[node.var_values] = node
            return node
        else:
            return self.nodes[node.var_values]

    def add_edge(self, edge):
        edge.src = self.add_node(edge.src)
        edge.dst = self.add_node(edge.dst)
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.edges:
            self.edges[edge_tuple] = edge
            self.nodes[edge.src.var_values].add_out_edge(edge)
            self.nodes[edge.dst.var_values].add_in_edge(edge)