from math import log
from z3 import Solver, Int, sat
from Search_Algos.BMC_GDFS.enabled import enabled
from Utils.misc import is_target, get_reaction_rate, get_total_outgoing_rate, trace_as_list
from Utils.encodings import state_encoding, bound_encoding, target_encoding, check_sat
from Graph.Edge import Edge
from Graph.Node import Node
import time

def BMC_GDFS(model, K, P, diag, target_var, target_index, target_value, E_count, sat_count, unsat_count, nodes):
    #initializing the unsat solver
    unsat_solver = Solver()
    unsat_solver.push()
    #
    
    #initializing the solver
    solver = Solver()
    solver.push()
    #

    #initial state of the model
    if model.get_initial_state() in nodes:
        init_node = nodes[model.get_initial_state()]
    else:
        init_node = Node()
        init_node.var_values = model.get_initial_state()
        init_node.initial_state = True
        nodes[init_node.var_values] = init_node
    #

    visited = set()
    stack = []
    stack.append((init_node, 0, enabled(model, diag, init_node, 0, P)))
    sat_time = 0
    unsat_time = 0

    while stack:
        n, p, E = stack[-1]
        
        if ((len(E)==0) or (is_target(n.var_values, target_index, target_value))):
            stack.pop()
            continue

        visited.add(n.var_values)
        r_index, to_sink, p_prime = E.pop()
        E_count = E_count + 1
        
        # n_prime <-- successor(n, r)
        dst_var_values = [None] * len(model.get_species_tuple())
        for j, coefficient in enumerate(model.get_reactions_vector()[r_index]):
            dst_var_values[j] = n.var_values[j] + coefficient
        dst_var_values = tuple(dst_var_values)
        if dst_var_values in nodes:
            n_prime = nodes[dst_var_values]
        else:
            n_prime = Node()
            n_prime.var_values = dst_var_values
            nodes[n_prime.var_values] = n_prime
        #

        if to_sink:
            rem_bound = (K - len(stack))
            if rem_bound in n_prime.unsat_bounds:
                continue
            
            check = True
            if rem_bound < 1:
                check = False
            
            time_temp = time.time()
            if check and check_sat(model, unsat_solver, n_prime.var_values, rem_bound, target_index, target_value):
                solver.pop()
                solver.push()
                solver.add(state_encoding(model, n_prime))
                for i in range (1, rem_bound + 1):
                    solver.add(bound_encoding(model, i))
                solver.add(target_encoding(rem_bound, target_var, target_value))

                if  solver.check()==sat:
                    #
                    sat_time = sat_time + (time.time() - time_temp)
                    sat_count = sat_count + 1
                    #

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
                    pi_list = trace_as_list(pi, rem_bound, model, diag.nodes)
                    flag_temp = True
                    for i, e in enumerate(pi_list):
                        diag.add_edge(e)
                        if e.src.var_values not in nodes:
                            nodes[e.src.var_values] = e.src
                        if e.dst.var_values not in nodes:
                            nodes[e.dst.var_values] = e.dst
                    #
                    
                    #add n_prime to stack
                    stack.append((n_prime, p_prime, enabled(model, diag, n_prime, p_prime, P)))
                    #
                    
                    for i in range(0, rem_bound-1):
                        if p_prime < P:
                            break
                        
                        if pi_list[i].dst.var_values in visited:
                            break
                        
                        node = pi_list[i].src
                        transition_probability = get_reaction_rate(node.var_values, model, pi_list[i].reaction)
                        transition_probability = float(transition_probability /get_total_outgoing_rate(node.var_values, model))
                        transition_probability = log(transition_probability, 10)
                        p_prime = p_prime + transition_probability
                        node = pi_list[i].dst
                        stack.append((node, p_prime, enabled(model, diag, node, p_prime, P)))
                else:
                    #
                    unsat_time = unsat_time + (time.time() - time_temp)
                    unsat_count = unsat_count + 1
                    #
                    n_prime.unsat_bounds.add(rem_bound)
                    
            else:
                if check:
                    #
                    unsat_time = unsat_time + (time.time() - time_temp)
                    unsat_count = unsat_count + 1
                    #
                    n_prime.unsat_bounds.add(rem_bound)

        else:
            if (n_prime.var_values not in visited):
                stack.append((n_prime, p_prime, enabled(model, diag, n_prime, p_prime, P)))
    
    return (diag, E_count, sat_count, unsat_count, sat_time, unsat_time)