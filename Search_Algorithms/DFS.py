import itertools
from Parser import Parser
from z3 import Solver, Int, And, Or, sat, simplify
from Graph import Node, Edge, Graph, check_probability
from Utils import get_total_outgoing_rate, get_reaction_rate, is_target
from math import log, pow, floor, ceil
import time
import random
import sys
import json

def CEX_GEN(json_data):
    start_time = time.time()

    #parsing parameters from json elements
    model_path = json_data['model_path']
    model_name = json_data['model_name']
    K = int(json_data['starting_bound'])
    csl_prop_lb = json_data['csl_property']
    prism_bin = json_data['prism_binary']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])

    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    print('Starting [DFS]... \n')

    diag = Graph()
  
    while (True):
        # file = open("Res/yeast/yeast_dfs.txt", "a")
        diag = DFS(model, diag, K)
        print('K = ' + str(K))
        print('\n')
        # file.write('K = ' + str(K))
        # file.write('\n')
        diag_size = len(diag.edges) + len(diag.nodes)
        print('number of transitions = ' + str(len(diag.edges)))
        print('\n')
        print('number of states = ' + str(len(diag.nodes)))
        print('\n')
        print('time stamp before model checking = ' + str(time.time()-start_time))
        print('\n')
        print('--')
        print('\n')
        # file.write('number of transitions = ' + str(len(diag.edges)))
        # file.write('\n')
        # file.write('number of states = ' + str(len(diag.nodes)))
        # file.write('\n')
        # file.write('time stamp before model checking = ' + str(time.time()-start_time))
        # file.write('\n')
        # file.write('--')
        # file.write('\n')
        result = check_probability(diag, model, model_name, prism_bin, csl_prop_lb)
        print('probability= ' + str(result))
        print('\n')
        print('size= ' + str(diag_size))
        print('\n')
        print('time= ' + str(time.time()-start_time))
        print('\n')
        # file.write('probability= ' + str(result))
        # file.write('\n')
        # file.write('size= ' + str(diag_size))
        # file.write('\n')
        # file.write('time= ' + str(time.time()-start_time))
        # file.write('\n')
        #print('='*50)
        K = K + 1
        diag = Graph()
        print("="*50)
        print('\n')
        # file.write("="*50)
        # file.write('\n')
        # file.close()

def DFS(model, diag, K):
    flag, init_node = diag.get_node(model.get_initial_state())
    if not flag:
        init_node = Node()
        init_node.var_values = model.get_initial_state()
        init_node.initial_state = True
    
    visited = {}
    stack = []
    stack.append((init_node, enabled(model, init_node), 0))

    while stack:
        n, E, depth = stack[-1]    
        # if ((len(E)==0) or (is_target(n.var_values, target_index, target_value))):
        if len(E) == 0 or depth > K:
            stack.pop()      
            continue

        visited[n.var_values] = depth
        r_index = E.pop()
        
        #
        # r_count[r_index] = r_count[r_index] + 1
        #

        # n_prime <-- successor(n, r)
        dst_var_values = [None] * len(model.get_species_tuple())
        for j, coefficient in enumerate(model.get_reactions_vector()[r_index]):
            dst_var_values[j] = n.var_values[j] + coefficient
        dst_var_values = tuple(dst_var_values)
        
        flag, n_prime = diag.get_node(dst_var_values)
        if not flag:
            n_prime = Node()
            n_prime.var_values = dst_var_values
            diag.add_node(n_prime)

        #add n->n_prime to diag
        edge_tuple = n.var_values + n_prime.var_values + (r_index,)
        if not diag.get_edge(edge_tuple)[0]:
            edge = Edge()
            edge.src = n
            edge.dst = n_prime
            edge.reaction = r_index
            edge.rate = get_reaction_rate(n.var_values, model, r_index)
            diag.add_edge(edge)
        #

        if (n_prime.var_values not in visited) or (depth + 1 < visited[n_prime.var_values]):
            stack.append((n_prime, enabled(model, n_prime), depth + 1))
    
    return diag

def enabled(model, curr_node):
    R = []
    for i, r in enumerate(model.get_reactions_vector()):
        rate = get_reaction_rate(curr_node.var_values, model, i) 
        if rate > 0:
           R.append(i)
    return R
