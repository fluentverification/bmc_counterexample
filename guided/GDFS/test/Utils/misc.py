from Graph.Node import Node
from Graph.Edge import Edge

def get_total_outgoing_rate(var_values, model):
    rate = 0
    for i in range(len(model.get_reactions_vector())):
        rate = rate + model.get_reaction_rate(var_values, i)
    return rate

def get_reaction_rate(var_values, model, reaction):
    return model.get_reaction_rate(var_values, reaction)

def is_target(var_values, target_index, target_value):
    if var_values[target_index] == target_value:
        return True
    return False

#turn a trace returned by z3 into an ordered list of edges [e0, e1, ...]
def trace_as_list(trace, bound, model, nodes):
    #a path is an ordered list of states(var assignments)
    ordered_variable_list = [None] * (bound+1)

    for d in trace.decls():
        index = int(d.name()[d.name().rfind('.')+1:])
        if index>bound:
            continue
        if ordered_variable_list[index] == None: 
            ordered_variable_list[index] = []
        ordered_variable_list[index].append([d.name()[:d.name().find('.')], int(str(trace[d]))])

    edge_list = [None] * bound

    for i in range(bound):
        src_var_values = [None] * len(model.get_species_tuple())
        dst_var_values = [None] * len(model.get_species_tuple())

        p1 = ordered_variable_list[i]
        p2 = ordered_variable_list[i+1]

        for e in p1:
            if 'selected_reaction' in e[0]:
                reaction = e[1]
                continue
            index = model.get_species_tuple().index(e[0])
            src_var_values[index] = e[1]
        
        for e in p2:
            if 'selected_reaction' in e[0]:
                continue
            index = model.get_species_tuple().index(e[0])
            dst_var_values[index] = e[1]

        if tuple(src_var_values) in nodes:
            src = nodes[tuple(src_var_values)]
        else:
            src = Node()
            src.var_values = tuple(src_var_values)
            if src.var_values == model.get_initial_state():
                src.initial_state = True
        
        if tuple(dst_var_values) in nodes:
            dst = nodes[tuple(dst_var_values)]
        else:
            dst = Node()
            dst.var_values = tuple(dst_var_values)

        edge = Edge()
        edge.src = src
        edge.dst = dst
        edge.reaction = reaction
        edge.rate = get_reaction_rate(src_var_values, model, reaction)
        edge_list[i] = edge
    
    return edge_list