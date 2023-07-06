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
    csl_prop_ub = json_data['csl_property_upper_bound']
    prism_bin = json_data['prism_binary']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    mc_step = int(json_data['model_check_step'])
    max_comb = int(json_data['max_combination'])

    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    print('Starting [Bounded DFS]... \n')
    #generate min_max dict for bound K
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)

    subsets = get_subsets(index_tuple, 1, max_comb)
    diag = Graph()

    while (True):
        min_max_dict = get_min_max(model, K, target_index, target_value, subsets) 
        diag = Bounded_DFS(model, min_max_dict, diag, target_var, target_index, target_value)
        print('K = ' + str(K))
        diag_size = len(diag.edges) + len(diag.nodes)
        print('diag_size = ' + str(diag_size))
        print('time stamp before model checking = ' + str(time.time()-start_time))
        print('----')
        result = check_probability(diag, model, model_name, prism_bin, csl_prop_lb)
        print('probability= ' + str(result))
        # result = check_probability(diag, model, model_name, prism_bin, csl_prop_ub)
        # print('upper-bound probability= ' + str(result))
        print('size= ' + str(diag_size))
        print('time= ' + str(time.time()-start_time))
        #print('='*50)
        K = K + 1

        dict_path = "/home/ubu/projects/probmc/bmc_counterexample/crn_verification/freq_dict.json"
        # r_dict = {'0' : 14.5, '1': 15}
        # with open(dict_path, 'w') as file:
        #     json.dump(r_dict, file)

        # if True and K == 15:
        #     freq_dict = modified_SSA(diag, model)
        #     # print(freq_dict)
        #     sum = 0
        #     for i, e in enumerate(freq_dict):
        #         sum = sum + freq_dict[i]
        #     for i, e in enumerate(freq_dict):
        #         freq_dict[i] = freq_dict[i] / sum
        #     print(freq_dict)
        #     print('##########')
        # k_shortest_paths(diag, model)
        # if K > 13 or False:
        #     astar(diag, model, model_name, prism_bin, csl_prop_lb)
        # if K == 15:
        #     freq_dict = modified_SSA(diag, model)
        #     print(freq_dict)
        #     print('##########')
        #     with open(dict_path, 'w') as file:
        #         json.dump(freq_dict, file)
        # if K > 15:
        #     break

        if K > 14:
            print(dijkstra(diag, model))
        diag = Graph()
        print ('=' * 50)
        
        



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


#utility functions

#returns all the subsets of a tuple with cardinality in range [L, U]
#returns a list of tuples where each tuple is a subset
def get_subsets (tuple, L, U):
    if U>len(tuple):
        U = len(tuple)
    return list(itertools.chain.from_iterable(itertools.combinations(tuple, r) for r in range(L, U+1)))

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
####################### SSA ##########################
######################################################

def modified_SSA(graph, model):
    returned_freq_dict = {}
    count = 0
    for i,r in enumerate(model.get_reactions_vector()):
            returned_freq_dict[i] = 0

    while (True):
        curr_state = graph.get_node(model.get_initial_state())[1]
        freq_dict = {}
        for i,r in enumerate(model.get_reactions_vector()):
            freq_dict[i] = 0
        time = 0
        flag = False
        while(True):
            #
            out_rate = 0
            for e in curr_state.out_edges.values():
                out_rate = out_rate + e.rate
            total_rate = out_rate
            #

            rand1 = random.uniform(0, 1)
            tau = float(1.0 / float(total_rate)) * log(float(1.0 / rand1))
            time = time + tau
            if time > 1000:
                break

            
            rand2 = random.uniform(0, out_rate)
            j = 0
            temp_sum = 0

            out_edges = [e for e in curr_state.out_edges.values()]
            while temp_sum <= rand2:
                temp_sum = temp_sum + out_edges[j].rate
                j = (j + 1)
            
            j = j - 1
            
            freq_dict[out_edges[j].reaction] = freq_dict[out_edges[j].reaction] + 1

            for i, e in enumerate(curr_state.out_edges.values()):
                if e.reaction == out_edges[j].reaction:
                    selected_edge = e
            curr_state = selected_edge.dst
            if curr_state.var_values[4] == 70:
                flag = True
                break
    
        # count = count + 1
        # for i, e in enumerate(freq_dict):
        #     returned_freq_dict[i] = returned_freq_dict[i] + freq_dict[i]
        # if count % 10 == 0:
        #     print(count)
        # if count > 1000:
        #     break
        if flag:
            count = count + 1
            sum = 0
            for i, e in enumerate(freq_dict):
                sum = sum + freq_dict[i]
            for i, e in enumerate(freq_dict):
                freq_dict[i] = freq_dict[i] / sum
            
            for i, e in enumerate(freq_dict):
                returned_freq_dict[i] = returned_freq_dict[i] + freq_dict[i]/1000
            print(count)
            if count >= 1000:
                break
    return returned_freq_dict



######################################################
######### error-trace generation Algorithms ##########
######################################################
import heapq

def astar(graph, model, model_name, prism_bin, csl_prop):
    max_rate = 0
    min_rate = sys.float_info.max
    for e in graph.edges:
        rate = 0 - log((graph.edges[e].rate) / (get_total_outgoing_rate(graph.edges[e].src.var_values, model)), 10)
        if rate > max_rate:
            max_rate = rate
        if rate < min_rate:
            min_rate = rate
    
    num_scores = 7
    score_list = []
    for i in range(num_scores+1):
        value = ((i)*(max_rate-min_rate)) / num_scores
        score_list.append(value)
    # quit()

    k = 10
    N = 1000
    
    #source of the search
    for n in graph.nodes:
        if graph.nodes[n].initial_state:
            source = graph.nodes[n]
    #
    
    queue = [(heuristic(source), 0, [source])]
    heapq.heapify(queue)
    paths = []

    while queue and len(paths) < k:
        _, cost, path = heapq.heappop(queue)
        current_node = path[-1]

        if current_node.var_values[4] == 70:
            paths.append((cost, path))

        else:
            neighbors = []
            for e in current_node.out_edges:
                neighbor = current_node.out_edges[e].dst
                total_rate = 0
                for e_2 in current_node.out_edges:
                    total_rate = total_rate + current_node.out_edges[e_2].rate
                weight = 0 - log((current_node.out_edges[e].rate) / (get_total_outgoing_rate(current_node.var_values, model)), 10)
                for i in range(len(score_list)-1):
                    if weight >= score_list[i] and weight <= score_list[i+1]:
                        weight = i + 1
                        break
                # print(weight)
                reaction = current_node.out_edges[e].reaction
            
                neighbors.append((neighbor, reaction, weight))
            
            for neighbor, reaction,  weight in neighbors:
                # if neighbor not in path[0::2]:
                #     new_cost = cost + weight
                #     new_path = path + [reaction] + [neighbor]
                #     heapq.heappush(queue, (new_cost + heuristic(neighbor), new_cost, new_path))
                new_cost = cost + weight
                new_path = path + [reaction] + [neighbor]
                heapq.heappush(queue, (new_cost + heuristic(neighbor), new_cost, new_path))
    # return 0
    print("^^^^^" * 5)
    print ("information on 1000* error-traces found by A* on the counterexample")
    ordered_paths = []
    for p in paths:
        #initializing reaction frequency and average population dictionaries
        reaction_dict = {}
        for i, r in enumerate(model.get_reactions_vector()):
            reaction_dict[i] = 0
        population_dict = {}
        for i, s in enumerate(model.get_species_tuple()):
            population_dict[s] = 0
        #
        path = p[1]
        path_graph = Graph()
        init_node = path[0]
        init_node.initial_state = True
        path_graph.add_node(init_node)
        for i, s in enumerate(population_dict):
            population_dict[s] = population_dict[s] + init_node.var_values[i]

        for i in range(1, len(path), 2):
            edge = Edge()
            edge.src = path[i-1]
            edge.dst = path[i+1]
            edge.reaction = path[i]
            edge.rate = get_reaction_rate(edge.src.var_values, model, edge.reaction)
            path_graph.add_edge(edge)
            reaction_dict[edge.reaction] = reaction_dict[edge.reaction] + 1
            for i, s in enumerate(population_dict):
                population_dict[s] = population_dict[s] + edge.dst.var_values[i]
        
        for i, e in enumerate(reaction_dict):
            reaction_dict[e] = reaction_dict[e] / floor(len(path)/ 2)
        for i, s in enumerate(population_dict):
            population_dict[s] = population_dict[s] / ceil(len(path) / 2)
        prob = check_probability(path_graph, model, model_name, prism_bin, csl_prop)
        # prob = 1
        print(prob)
        ordered_paths.append([prob, reaction_dict, population_dict])
    
    ordered_paths =  sorted(ordered_paths, key=lambda x: x[0], reverse=True)
    # ordered_paths = ordered_paths[:N]
    
    reaction_dict_wa = {}
    for i, r in enumerate(model.get_reactions_vector()):
        reaction_dict_wa[i] = 0
    population_dict_wa = {}
    for i, s in enumerate(model.get_species_tuple()):
        population_dict_wa[s] = 0
    
    total_weight = 0
    for i, p in enumerate(ordered_paths):
        weight = p[0] 
        curr_reaction_dict = p[1]
        curr_population_dict = p[2]
        for j, r in enumerate(curr_reaction_dict):
            reaction_dict_wa[r] = reaction_dict_wa[r] + (curr_reaction_dict[r] * weight)
        for j, s in enumerate(curr_population_dict):
            population_dict_wa[s] = population_dict_wa[s] + (curr_population_dict[s] * weight)
        total_weight = total_weight + weight

    for j, r in enumerate(reaction_dict_wa):
        reaction_dict_wa[r] = reaction_dict_wa[r] / total_weight
    for j, s in enumerate(population_dict_wa):
        population_dict_wa[s] = population_dict_wa[s] / total_weight
        
    print("number of error-traces : " + str(len(ordered_paths)))
    print("average population of spicies among error-traces: ")
    print(population_dict_wa)
    print("frequency of reactions among error-traces: ")
    print(reaction_dict_wa)
    print("^^^^^" * 5)

def heuristic(node):
    # Replace with your own heuristic function
    # return 0
    return floor((70 - node.var_values[4])/10)
######################################################
######################################################
######################################################


######################################################
######################Dijkstra########################
######################################################
import heapq

def dijkstra(graph, model):
    print(model.get_species_tuple())
    target_node = Node()
    target_node.var_values = (-2, -2, -2, -2, -2)
    t_edges = []
    for n in graph.nodes:
        if graph.nodes[n].var_values[1] == 0:
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
    current_node = (-2, -2, -2, -2, -2)
    path = []
    while current_node:
        path.insert(0, current_node)
        edge = previous_edges[current_node]
        if edge:
            path.insert(0, edge)
        current_node = previous_nodes[current_node]

    return path


# def shortest_path(graph, start_node, end_node):
#     distances, previous_nodes, previous_edges = dijkstra(graph, start_node)
#     path = []
#     current_node = end_node

#     # Trace back the path from end_node to start_node using the previous nodes and edges
#     while current_node:
#         path.insert(0, current_node)
#         edge = previous_edges[current_node]
#         if edge:
#             path.insert(0, edge)
#         current_node = previous_nodes[current_node]

#     return path

# graph = {
#     'A': {'B': 4, 'C': 2, 'C': 1},
#     'B': {'D': 5, 'E': 2},
#     'C': {'B': 1, 'D': 8},
#     'D': {'E': 3},
#     'E': {}
# }

# start_node = 'A'
# end_node = 'E'

# path = shortest_path(graph, start_node, end_node)
# print("Shortest path from", start_node, "to", end_node, ":", path)


# quit()