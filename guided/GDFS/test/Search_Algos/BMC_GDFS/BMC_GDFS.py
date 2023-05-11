from math import log
from z3 import Solver, Int, sat
from Search_Algos.BMC_GDFS.enabled import enabled
from Utils.misc import is_target, get_reaction_rate, get_total_outgoing_rate, trace_as_list
from Utils.encodings import initial_state_encoding, bound_encoding, target_encoding
from Graph.Edge import Edge
from Graph.Node import Node

def BMC_GDFS(model, K, P, diag, target_var, target_index, target_value):
    #initializing the solver
    solver = Solver()
    solver.add(initial_state_encoding(model))
    for i in range(1, K+1):
        solver.add(bound_encoding(model, i))
    solver.add(target_encoding(K, target_var, target_value))
    solver.push()
    #

    #initial state of the model
    flag, node = diag.get_node(model.get_initial_state())
    if flag:
        init_node = node
    else:
        init_node = Node()
        init_node.var_values = model.get_initial_state()
        init_node.initial_state = True
    #

    visited = set()
    stack = []
    stack.append((init_node, 0, enabled(model, diag, init_node, 0, P)))

    while stack:
        n, p, E = stack[-1]
        
        if ((len(E)==0) or (False) 
            or (is_target(n.var_values, target_index, target_value))):
            stack.pop()
            solver.pop()
            continue

        visited.add(n.var_values)
        r_index, to_sink, p_prime = E.pop()
        
        # n_prime <-- successor(n, r)
        dst_var_values = [None] * len(model.get_species_tuple())
        for j, coefficient in enumerate(model.get_reactions_vector()[r_index]):
            dst_var_values[j] = n.var_values[j] + coefficient
        dst_var_values = tuple(dst_var_values)
        if diag.get_node(dst_var_values)[0]:
            n_prime = diag.get_node(dst_var_values)[1]
        else:
            n_prime = Node()
            n_prime.var_values = dst_var_values
        #

        if to_sink:
            x = Int('selected_reaction.'+ str(len(stack)-1))
            if solver.check((x==r_index))==sat:
                #add n->n_prime to diag
                edge = Edge()
                edge.src = n
                edge.dst = n_prime
                edge.reaction = r_index
                edge.rate = get_reaction_rate(n.var_values, model, r_index)
                diag.add_edge(edge)
                #

                #add pi to graph
                pi = solver.model()
                pi_list = trace_as_list(pi, K, model, diag.nodes)
                for i, e in enumerate(pi_list):
                    diag.add_edge(e)
                    if e.src.var_values == n_prime.var_values:
                        index = i
                #
                
                for i in range(index, K):
                    if p_prime < P:
                        break
                    
                    if pi_list[i].src.var_values in visited:
                        break
                    
                    x = Int('selected_reaction.'+ str(len(stack)-1))
                    solver.add(x == r_index)
                    solver.push()
                    
                    node = pi_list[i].src
                    transition_probability = get_reaction_rate(node.var_values, model, pi_list[i].reaction)
                    transition_probability = float(transition_probability /get_total_outgoing_rate(node.var_values, model))
                    transition_probability = log(transition_probability, 10)
                    p_prime = p_prime + transition_probability
                    stack.append((node, p_prime, enabled(model, diag, node, p_prime, P)))
                    
                    

        else:
            if n_prime.var_values not in visited:
                x = Int('selected_reaction.'+ str(len(stack)-1))
                solver.add(x == r_index)
                solver.push()
                stack.append((n_prime, p_prime, enabled(model, diag, n_prime, p_prime, P)))
    
    return diag