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
    mc_step = int(json_data['model_check_step'])
    max_comb = int(json_data['max_combination'])
    max_comb = 1

    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    print('Starting [Bounded DFS]... \n')
    #generate min_max dict for bound K
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)

    subsets = get_subsets(index_tuple, 1, max_comb)
    # subsets.append(index_tuple)

    diag = Graph()

    mc_time = 0
    construction_time = 0
    constraint_time = 0    
    while (True):
        file = open("./Res/circuit/circuit_single_species_2.txt", "a")
        stamp = time.time()
        min_max_dict = get_min_max(model, K, target_index, target_value, subsets)
        constraint_time = constraint_time + (time.time() - stamp)
        file.write("time stamp after constraint generation: " + str(time.time()-start_time))
        file.write("\n")
        stamp = time.time()
        diag = Bounded_DFS(model, min_max_dict, diag, target_var, target_index, target_value)
        construction_time = construction_time + (time.time() - stamp)
        file.write('K = ' + str(K))
        file.write("\n")
        diag_size = len(diag.edges) + len(diag.nodes)
        file.write("number of transitions: " + str(len(diag.edges)))
        file.write("\n")
        file.write("number of states: " + str(len(diag.nodes)))
        file.write("\n")
        # print('diag_size = ' + str(diag_size))
        file.write('time stamp before model checking = ' + str(time.time()-start_time))
        file.write("\n")
        file.write('----')
        file.write("\n")
        stamp = time.time()
        result = check_probability(diag, model, model_name, prism_bin, csl_prop_lb)
        mc_time = mc_time + (time.time() - stamp)
        file.write('probability= ' + str(result))
        file.write("\n")
        file.write('size= ' + str(diag_size))
        file.write("\n")
        file.write('time= ' + str(time.time()-start_time))
        file.write("\n")
        file.write("total constraint generation time= " + str(constraint_time))
        file.write("\n")
        file.write("total graph construction time= " + str(construction_time))
        file.write("\n")
        file.write("total model checking time= " + str(mc_time))
        file.write("\n")
        file.write("="*50)
        file.write("\n")
        file.close()
        #print('='*50)
        K = K + 1

        # print('time stamp before generating traces ' + str(time.time()-start_time))
        # print(dijkstra(diag, model, target_index, target_value))
        # print('time stamp after generating traces ' + str(time.time()-start_time))


        diag = Graph()
        # print('=' * 50)
        
        



def Bounded_DFS(model, min_max_dict, diag, target_var, target_index, target_value):
    # r_count = {}
    # for i, r in enumerate(model.get_reactions_vector()):
    #     r_count[i] = 0

    flag, init_node = diag.get_node(model.get_initial_state())
    if not flag:
        init_node = Node()
        init_node.var_values = model.get_initial_state()
        init_node.initial_state = True
    
    visited = set()
    stack = []
    stack.append((init_node, enabled(model, init_node)))

    while stack:
        n, E = stack[-1]    
        if ((len(E)==0) or (is_target(n.var_values, target_index, target_value))):
        # if len(E) == 0:
            stack.pop()      
            continue

        visited.add(n.var_values)
        r_index = E.pop()
        
        #
        # r_count[r_index] = r_count[r_index] + 1
        #

        # n_prime <-- successor(n, r)
        dst_var_values = [None] * len(model.get_species_tuple())
        for j, coefficient in enumerate(model.get_reactions_vector()[r_index]):
            dst_var_values[j] = n.var_values[j] + coefficient
        dst_var_values = tuple(dst_var_values)
        
        if check_state(dst_var_values, min_max_dict):
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

            if n_prime.var_values not in visited:
                stack.append((n_prime, enabled(model, n_prime)))
    
    #
    # sum = 0
    # for v in r_count.values():
    #     sum = sum + v
    # r_freq = {}
    # for i, r in enumerate(r_count):
    #     r_freq[r] = float(r_count[r] / sum)
    # return diag, r_count, r_freq
    
    
    return diag

def DFS(model, min_max_dict, diag, target_var, target_index, target_value, K):
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
        if len(E) == 0 or depth > K or (is_target(n.var_values, target_index, target_value)):
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

        if check_state(dst_var_values, min_max_dict):
            if (n_prime.var_values not in visited) or (depth + 1 < visited[n_prime.var_values]):
                stack.append((n_prime, enabled(model, n_prime), depth + 1))
    
    return diag

#utility functions

############ Functions Generating Min_Max via SMT-Solving ############
#returns all the subsets of a tuple with cardinality in range [L, U]
#returns a list of tuples where each tuple is a subset
def get_subsets (tuple, L, U):
    if U>len(tuple):
        U = len(tuple)
    return list(itertools.chain.from_iterable(itertools.combinations(tuple, r) for r in range(L, U+1)))
#

#generating a dictionary. 
#Keys: combinations of variables:(s1,), (s1,s2), ...
#Values: a list with l[0] the minimum value and l[1] the maximum value for that key
def get_min_max(model, bound, target_index, target_value, subsets):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in subsets:
        min_max[s] = [-1, -1]
    
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_i_" + str(i))
        vars.append(x)
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_t_" + str(i))
        vars.append(x)
    constraints = []
    
    #first constraint (n_i_0,n_t_0,n_i_1,n_t_1,... >=0)
    for i in vars:
        constraints.append(i>=0)
    
    #second constraint (n_i_0+n_i_t+n_i_1+...<=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum<=bound))
    
    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            vars.append([Int("n_i_" + str(i)), r[target_index]])
            vars.append([Int("n_t_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    #fourth constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n_i_"+str(j)), r[i]])
                vars.append([Int("n_t_"+str(j)), r[i]])
                vars_i.append([Int("n_i_"+str(j)), r[i]])
        sum = model.get_initial_state()[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)
        sum = model.get_initial_state()[i]
        for j in vars_i:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)

    solver = Solver()
    solver.add(And(constraints))

    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            max_value = 0
            for j in e:
                max_value = max_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    max_value = max_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum>max_value)
        solver.pop()
    
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            min_value = 0
            for j in e:
                min_value = min_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    min_value = min_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][0] = min_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum<min_value)
        solver.pop()
    return min_max
#

def check_state(var_values, min_max):
    for e in min_max:
        value = 0
        for j in e:
            value = value + var_values[j]
        if value > min_max[e][1] or value < min_max[e][0]:
            return False
    return True


def enabled(model, curr_node):
    R = []
    for i, r in enumerate(model.get_reactions_vector()):
        rate = get_reaction_rate(curr_node.var_values, model, i) 
        if rate > 0:
           R.append(i)
    return R


######################################################
######################Dijkstra########################
######################################################
import heapq

def dijkstra(graph, model, target_index, target_value):
    print(model.get_species_tuple())
    target_node = Node()
    target_node.var_values = (-2, -2, -2, -2, -2, -2, -2)
    t_edges = []
   
   
    
    
   
    for n in graph.nodes:
        
        AmtR = graph.nodes[n].var_values[0]
        BetI = graph.nodes[n].var_values[1]
        HlyIIR = graph.nodes[n].var_values[2]
        PhIF = graph.nodes[n].var_values[3]
        YFP = graph.nodes[n].var_values[4]

        # if is_target(graph.nodes[n].var_values, target_index, target_value):
        if  AmtR >=220 and BetI>=70 and PhIF >= 70 and YFP >= 60 and HlyIIR>=280:
            edge = Edge()
            edge.reaction = -2
            edge.src = graph.nodes[n]
            edge.dst = target_node
            edge.rate = -2
            #graph.add_edge(edge)
            t_edges.append(edge)
    for edge in t_edges:
        graph.add_edge(edge)

    for n in graph.nodes:
        if graph.nodes[n].initial_state:
            start_node = n
    # Initialize distances, previous nodes, and previous edges dictionaries
    distances = {node: float('inf') for node in graph.nodes}
    distances[start_node] = 0
    previous_nodes = {node: None for node in graph.nodes}
    previous_edges = {node: None for node in graph.nodes}

    # Priority queue to store nodes with their respective distances
    priority_queue = [(0, start_node)]

    while priority_queue:
        current_distance, current_node = heapq.heappop(priority_queue)

        # If the current distance is greater than the known distance, skip
        if current_distance > distances[current_node]:
            continue

        # Visit all neighboring nodes and update their distances
        neighbors = []
        for e in graph.nodes[current_node].out_edges:
            edge = graph.edges[e]
            neighbor = edge.dst.var_values
            reaction = edge.reaction
            if edge.reaction == -2:
                weight = 0
            else:
                total_rate = get_total_outgoing_rate(edge.src.var_values, model)
                weight = 0 - log(edge.rate / total_rate ,10)
            neighbors.append((neighbor, weight, reaction))
        for neighbor, weight, reaction in neighbors:
            distance = current_distance + weight

            # If a shorter distance is found, update the distance, previous node, and previous edge
            if distance < distances[neighbor]:
                distances[neighbor] = distance
                previous_nodes[neighbor] = current_node
                previous_edges[neighbor] = reaction
                heapq.heappush(priority_queue, (distance, neighbor))

    #return distances, previous_nodes, previous_edges
    current_node = (-2, -2, -2, -2, -2, -2, -2)
    path = []
    while current_node:
        path.insert(0, current_node)
        edge = previous_edges[current_node]
        if edge:
            path.insert(0, edge)
        current_node = previous_nodes[current_node]

    return path