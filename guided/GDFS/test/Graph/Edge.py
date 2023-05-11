class Edge:
    def __init__(self):
        self.src = None
        self.dst = None
        self.rate = 0
        self.reaction = -1

    def __eq__(self, edge) -> bool:
        src_check = (self.src.var_values == edge.src.var_values)
        dst_check = (self.dst.var_values == edge.dst.var_values)
        reaction_check = (self.reaction == edge.reaction)
        return src_check and dst_check and reaction_check
    
    def get_tuple(self):
        return self.src.var_values + self.dst.var_values + (self.reaction,)
