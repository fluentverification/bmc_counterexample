from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, get_initial_state, get_encoding, backtrack
from Graph import Node, Edge, Graph
import time
from z3 import Int, And, Solver, sat, Or
from math import log


###################################################################
#turn a trace returned by z3 into an ordered list of edges
def trace_as_list(trace, bound, model):
    #a path is an ordered list of states(var assignments)
    ordered_variable_list = [None] * (bound+1)

    for d in trace.decls():
        index = int(d.name()[d.name().rfind('.')+1:])
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

        src = Node()
        src.var_values = src_var_values
        dst = Node()
        dst.var_values = dst_var_values

        edge = Edge()
        edge.src = src
        edge.dst = dst
        edge.reaction = reaction
        edge.rate = get_reaction_rate(src_var_values, model, reaction)
        edge_list[i] = edge
    
    return edge_list
###################################################################

###### Greedy DFS guided by BMC and a probability threshold #######

def BMC_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, threshold, start_time, diag):  
    nodes = {}
    diag_size = len(diag.nodes) + len(diag.edges)

    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()

    ### Setting up always necessary constraints for to found an error trace with length = bound ###
    solver = Solver()
    solver.add(get_initial_state(model))

    for j in range(1, bound+1):
        solver.add(get_encoding(model, j))

    x = Int(target_var + '.' + str(bound))
    property_constraint = (x==target_value)
    solver.add(property_constraint)
    solver.push()
    ################################################################################################

    visited = set()
    stack = []
    succ = []

    for i, r in enumerate(model.get_reactions_vector()):
        rate = get_reaction_rate(init_node.var_values, model, i)
        total_rate = get_total_outgoing_rate(init_node.var_values, model)
        if rate > 0:
            dst_var_values = [None] * len(model.get_species_tuple())
            for j, coefficient in enumerate(r):
                dst_var_values[j] = init_node.var_values[j] + coefficient
            dst_var_values = tuple(dst_var_values)

            if dst_var_values in nodes:
                succ_node = nodes[dst_var_values]
            else:
                succ_node = Node()
                succ_node.var_values = dst_var_values
                nodes[succ_node.var_values] = succ_node
            succ.append((succ_node, log(rate/total_rate, i)))
    
    succ = sorted(succ, key=lambda x: x[1])
    stack.push((init_node, 0, succ, 0))
    current_depth = 0

    while stack:
        current_state, path_probability, succ, depth = stack[-1]
        if current_state.var_values in visited:
            stack.pop()
            continue      
        visited.add(current_state.var_values)
        
        if len(succ) == 0:
            stack.pop()
            continue
        
        if (is_target(current_state.var_values, target_index, target_value) 
            or diag.get_node(current_state.var_values)[0]):
            if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                print('probability= ' + str(result))
                diag_size = len(diag.edges) + len(diag.nodes)
                print('size= ' + str(diag_size))
                print('time= ' + str(time.time()-start_time))
                print('='*30)
                if time.time()-start_time > 3600:
                    quit()

        n, p, r = succ.pop()

        if depth < current_depth:
            for i in range (current_depth-depth):
                solver.pop()
        
        x = Int('selected_reaction.' + str(depth - 1))
        solver.add(x == r)
        solver.push()

        if solver.check()==sat:
            pi = solver.model()
            pi_list = trace_as_list(pi, bound, model)
            for e in pi_list:
                diag.add_edge(e)
            
            for i in range(len(pi_list)):
                if current_state.var_values == pi_list[i].src.var_values:
                    index = i
                    break
            
            for i in range(index, bound):
                if i == index:
                    p_prime = p
                else:
                    total_outgoing_rate = get_total_outgoing_rate(current_state.var_values, model)
                    p_prime = p + log(pi_list[i].rate / total_outgoing_rate)
                if p_prime < threshold:
                    break
                
                current_state = pi_list[i].src
                succ = []
                
                for i, r in enumerate(model.get_reactions_vector()):
                    rate = get_reaction_rate(current_state.var_values, model, i)
                    total_rate = get_total_outgoing_rate(current_state.var_values, model)
                    if rate > 0:
                        dst_var_values = [None] * len(model.get_species_tuple())
                        for j, coefficient in enumerate(r):
                            dst_var_values[j] = current_state_node.var_values[j] + coefficient
                        dst_var_values = tuple(dst_var_values)

                        if dst_var_values in nodes:
                            succ_node = nodes[dst_var_values]
                        else:
                            succ_node = Node()
                            succ_node.var_values = dst_var_values
                            nodes[succ_node.var_values] = succ_node
                        succ.append((succ_node, log(rate/total_rate, i)))
                
                depth = depth + 1
                x = Int('selected_reaction.' + str(depth - 1))
                solver.add(x == pi_list[i].reaction)
                solver.push()
                stack.push(current_state, p_prime, succ, depth)
    return diag