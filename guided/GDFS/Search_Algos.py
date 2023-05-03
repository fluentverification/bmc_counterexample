from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, check_sat, get_score, get_initial_state, get_encoding
from Graph import Node, Edge, Graph
from heapq import heapify, heappush, heappop
import time
import numpy as np
from z3 import Int, And, Solver, sat, Or
from math import log

def construct_trace(graph, node):
    curr_node = node
    flag = True
    while((not graph.get_node(curr_node.var_values)[0]) or flag):
        flag = False
        graph.add_edge(curr_node.incoming_edge)
        curr_node = curr_node.incoming_edge.src

def post_process(graph, target_index, target_value):
    new_graph = Graph()
    visited = {}
    stack = []
    for n in graph.nodes.values():
        if is_target(n.var_values, target_index, target_value):
            new_graph.add_node(n)
            stack.append(n)
    while stack:
        curr_node = stack.pop()
        visited[curr_node.var_values] = curr_node
        for e in curr_node.in_edges.values():
            new_graph.add_edge(e)
            if not e.src.var_values in visited:
                stack.append(e.src)
    return new_graph

#guided greedy DFS where all the nodes allowed by the minmax dictionary are added to search stack. 
#All nodes added to the search stack are added to the graph.
def ggdfs(graph, start_node, model, target_index, target_value, min_max):
    flag = False
    visited = {}
    stack = [start_node]
    while stack:
        node = stack.pop()
        if not node.var_values in visited:
            visited[node.var_values] = node
            if is_target(node.var_values, target_index, target_value):
                flag = True
                neighbors = []
            else:
                neighbors = []
                for i, r in enumerate(model.get_reactions_vector()):
                    rate = get_reaction_rate(node.var_values, model, i)
                    total_rate = get_total_outgoing_rate(node.var_values, model)
                    if rate > 0:
                        dst_var_values = [None] * len(model.get_species_tuple())
                        for j, coefficient in enumerate(r):
                            dst_var_values[j] = node.var_values[j] + coefficient
                        dst_var_values = tuple(dst_var_values)
                        if check_sat(dst_var_values, min_max):
                            new_node = Node()
                            new_node.var_values = dst_var_values
                            new_node.reachability_probability = float(rate/total_rate)*node.reachability_probability
                            edge = Edge()
                            edge.src = node
                            edge.dst = new_node
                            edge.rate = rate
                            edge.reaction = i

                            if not graph.get_edge(edge.get_tuple())[0]:
                                graph.add_edge(edge)

                            neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append(neighbor)
    return flag
    
#guided greedy DFS where all the nodes allowed by the minmax dictionary are added to search stack. 
#All nodes added to the search stack are added to the graph.
# def ggdfs_prob(graph, start_node, model, target_index, target_value, min_max, prob_thresh, size):
#     flag = False
#     visited = {}
#     stack = [(start_node, 1.0)]
#     while stack:
#         node, prob = stack.pop()
        
#         if not node.var_values in visited:
#             visited[node.var_values] = node
#             if is_target(node.var_values, target_index, target_value):
#                 flag = True
#                 neighbors = []
#             else:
#                 neighbors = []
#                 for i, r in enumerate(model.get_reactions_vector()):
#                     rate = get_reaction_rate(node.var_values, model, i)
#                     total_rate = get_total_outgoing_rate(node.var_values, model)
#                     if rate > 0:
#                         dst_var_values = [None] * len(model.get_species_tuple())
#                         for j, coefficient in enumerate(r):
#                             dst_var_values[j] = node.var_values[j] + coefficient
#                         dst_var_values = tuple(dst_var_values)
#                         if check_sat(dst_var_values, min_max) and prob>prob_thresh:
#                             new_node = Node()
#                             new_node.var_values = dst_var_values
#                             new_node.reachability_probability = (float(rate/total_rate)*prob)
#                             edge = Edge()
#                             edge.src = node
#                             edge.dst = new_node
#                             edge.rate = rate
#                             edge.reaction = i

#                             if prob>prob_thresh:
#                                 if new_node.var_values in visited:
#                                     old_prob = visited[new_node.var_values].reachability_probability
#                                     if new_node.reachability_probability > old_prob:
#                                         del visited[new_node.var_values]

#                                 # if new_node.var_values in visited:
#                                 #     del visited[new_node.var_values]
#                                 if not graph.get_edge(edge.get_tuple())[0]:
#                                     graph.add_edge(edge)
#                                 neighbors.append(new_node)
#             neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
#             for neighbor in neighbors:
#                 if neighbor.var_values not in visited:
#                     stack.append((neighbor, new_node.reachability_probability))
#     return flag



def ggdfs_prob(graph, start_node, model, target_index, target_value, min_max, prob_thresh, size):
    flag = False
    visited = {}
    stack = [(start_node, 1.0)]

    while stack:
        node, prob = stack.pop()
        if node.var_values not in visited:
            visited[node.var_values] = node
            if is_target(node.var_values, target_index, target_value):
                flag = True
                neighbors = []
            else:
                neighbors = []
                for i, r in enumerate(model.get_reactions_vector()):
                    rate = get_reaction_rate(node.var_values, model, i)
                    total_rate = get_total_outgoing_rate(node.var_values, model)
                    if rate > 0 :
                        dst_var_values = [None] * len(model.get_species_tuple())
                        for j, coefficient in enumerate(r):
                            dst_var_values[j] = node.var_values[j] + coefficient
                        dst_var_values = tuple(dst_var_values)
                        if check_sat(dst_var_values, min_max) and ((prob*float(rate/total_rate))>=prob_thresh):
                            new_node = Node()
                            new_node.var_values = dst_var_values
                            edge = Edge()
                            edge.src = node
                            edge.dst = new_node
                            edge.rate = rate
                            edge.reaction = i

                            if not graph.get_edge(edge.get_tuple())[0]:
                                    graph.add_edge(edge)
                            neighbors.append((new_node, prob*float(rate/total_rate)))
                
                neighbors = sorted(neighbors, key=lambda n: n[1], reverse = True)
                for neighbor in neighbors:
                    stack.append(neighbor)
    return flag

######################################################
########### Best First Search Methods ################

def h(reach_prob, var_values, model, min_max_list):
    #return -np.float128(reach_prob*np.float128(get_score(var_values, min_max_list)))
    return -np.float128(reach_prob)

def backtrack(graph, node, model):
    visited = set()
    stack = [node]
    while (stack):
        curr_node = stack.pop()
        if curr_node.var_values not in visited:
            for e in curr_node.in_edges.values():
                if e.src.var_values not in graph.nodes:   
                    stack.append(e.src)
                graph.add_edge(e)

        visited.add(curr_node.var_values)

def backtrack_2(graph, node, model):
    visited = set()
    stack = [node]
    while (stack):
        curr_node = stack.pop()
        if curr_node.var_values not in visited:
            for e in curr_node.in_edges.values():
                if e.src.var_values != model.get_initial_state():   
                    stack.append(e.src)
                graph.add_edge(e)

        visited.add(curr_node.var_values)
            


def XBF(model, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, min_max_list):
    start_time = time.time()
    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()
    init_node.reachability_probability = np.float128(1.0)
    nodes = {}
    nodes[init_node.var_values] = init_node
    # diag <-- empty graph
    # closed <-- empty hash table
    # open_ <-- empty priority queue (implemented using heap functions)
    # map_ <-- mapping node's var_values to a list of size 2 containing the hueristic value and the node object
    # map_ gives access to the list elements of open_. after an element of the list of size 2 is changed, the 
    # open_ should get heapified again in order to maintain the logic of the heap
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)
    open_ = []
    closed = {}
    map_ = {}

    node_list = [h(1.0, init_node.var_values, model, min_max_list), init_node]
    map_[init_node.var_values] = node_list
    heappush(open_, node_list)

    while (open_):
        h_s = heappop(open_)
        old_h = h_s[0]
        s = h_s[1]
        del map_[s.var_values]
        closed[s.var_values] = [old_h, s]

        
        for i, r in enumerate(model.get_reactions_vector()):
            rate = get_reaction_rate(s.var_values, model, i)
            total_rate = get_total_outgoing_rate(s.var_values, model)
            
            #
            dst_var_values = [None] * len(model.get_species_tuple())
            for j, coefficient in enumerate(r):
                dst_var_values[j] = s.var_values[j] + coefficient
            #

            if rate > 0 and check_sat(dst_var_values, min_max_list):
                dst_var_values = [None] * len(model.get_species_tuple())
                for j, coefficient in enumerate(r):
                    dst_var_values[j] = s.var_values[j] + coefficient
                dst_var_values = tuple(dst_var_values)
                if dst_var_values in nodes:
                    s_prime = nodes[dst_var_values]
                else:
                    s_prime = Node()
                    s_prime.var_values = dst_var_values
                    nodes[s_prime.var_values] = s_prime
                edge = Edge()
                edge.src = s
                edge.dst = s_prime
                edge.rate = rate
                edge.reaction = i
                s_prime.add_in_edge(edge)
                if (np.multiply(s.reachability_probability, np.float128(rate/total_rate))) > s_prime.reachability_probability:
                    s_prime.reachability_probability = np.multiply(s.reachability_probability, np.float128(rate/total_rate))

                new_h = h(s_prime.reachability_probability, s_prime.var_values, model, min_max_list)

                if (is_target(s_prime.var_values, target_index, target_value)
                    or diag.get_node(s_prime.var_values)[0]):
                    backtrack(diag, s_prime, model)
                    if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                        result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                        print('probability= ' + str(result))
                        diag_size = len(diag.edges) + len(diag.nodes)
                        print('size= ' + str(diag_size))
                        print('time= ' + str(time.time()-start_time))
                        print('='*50)
                        if time.time()-start_time > 360000:
                            quit()
            
                if ((s_prime.var_values not in map_) 
                    and (s_prime.var_values not in closed)
                    and (not is_target(s_prime.var_values, target_index, target_value))):
                    node_list = [new_h, s_prime]
                    heappush(open_, node_list)
                    map_[s_prime.var_values] = node_list
                else:
                    if (s_prime.var_values in map_):
                        if new_h < map_[s_prime.var_values][0]:
                            map_[s_prime.var_values][0] = new_h
                            #heapify(open_)
                    elif (s_prime.var_values in closed):
                        if new_h < closed[s_prime.var_values][0]:
                            del closed[s_prime.var_values]
                            node_list = [new_h, s_prime]
                            heappush(open_, node_list)
                            if s_prime.var_values in map_:
                                raise Exception("a node is both in open and closed")
                            map_[s_prime.var_values] = node_list

######################################################
########### Greedy DFS guided by BMC #################

def Guided_Greedy_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound):
    start_time = time.time()
    bound = init_bound
    old_bound = bound
    nodes = {}
    
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)

    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()
    solver = Solver()
    solver.add(get_initial_state(model))

    for j in range(1, bound+1):
        solver.add(get_encoding(model, j))

    solver.push()
    x = Int(target_var + '.' + str(bound))
    property_constraint = (x==target_value)
    solver.add(property_constraint)  

    while True:
        current_depth = 0
        print('+++++' * 10)
        print('bound = ' + str(bound) + ' \n')
        if bound == (old_bound+1):
            solver.pop()
            solver.add(get_encoding(model, bound))
            solver.push()
            x = Int(target_var + '.' + str(bound))
            property_constraint = (x==target_value)
            solver.add(property_constraint)
            old_bound = bound        
        elif (bound != init_bound):
            raise Exception("bound did not increase by one after an iteration")

        visited = {}
        stack = [[init_node, 1.0, 0]]
        last_visited_depth = 0
        
        while stack:
            node, reachability_probability, depth = stack.pop()

            if depth < current_depth:
                for i in range(0, current_depth - depth):
                    solver.pop()
                current_depth = depth
                last_visited_depth = depth
            
            if depth != 0:
                if depth != last_visited_depth:
                    solver.push()
                    x = Int('selected_reaction.' + str(depth - 1))
                    solver.add(x == node.last_in_edge.reaction)
                    last_visited_depth = depth
                    
                else:
                    solver.pop()
                    solver.push()
                    x = Int('selected_reaction.' + str(depth - 1))
                    solver.add(x == node.last_in_edge.reaction)
            
            if depth > current_depth:
                raise Exception("current_depth is behind the depth of the node at the top of the stack")

            if node.var_values not in visited:
                visited[node.var_values] = node
                if is_target(node.var_values, target_index, target_value):
                    neighbors = []
                    backtrack(diag, node, model)
                    if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                            result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                            print('bound= ' + str(bound))
                            print('probability= ' + str(result))
                            diag_size = len(diag.edges) + len(diag.nodes)
                            print('size= ' + str(diag_size))
                            print('time= ' + str(time.time()-start_time))
                            print('='*50)
                            if time.time()-start_time > 3600:
                                quit()
                else:
                    if (diag.get_node(node.var_values)[0]):
                        backtrack(diag, node, model)
                        if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                            result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                            print('bound= ' + str(bound))
                            print('probability= ' + str(result))
                            diag_size = len(diag.edges) + len(diag.nodes)
                            print('size= ' + str(diag_size))
                            print('time= ' + str(time.time()-start_time))
                            print('='*50)
                            if time.time()-start_time > 3600:
                                quit()
                    
                    neighbors = []
                    if solver.check()==sat:
                        for i, r in enumerate(model.get_reactions_vector()):
                            rate = get_reaction_rate(node.var_values, model, i)
                            total_rate = get_total_outgoing_rate(node.var_values, model)
                            if rate > 0:
                                dst_var_values = [None] * len(model.get_species_tuple())
                                for j, coefficient in enumerate(r):
                                    dst_var_values[j] = node.var_values[j] + coefficient
                                dst_var_values = tuple(dst_var_values)

                                if dst_var_values in nodes:
                                    new_node = nodes[dst_var_values]
                                else:
                                    new_node = Node()
                                    new_node.var_values = dst_var_values
                                    nodes[new_node.var_values] = new_node
                                    
                                
                                
                                edge = Edge()
                                edge.src = node
                                edge.dst = new_node
                                edge.rate = rate
                                edge.reaction = i
                                
                                new_node.add_in_edge(edge)
                                new_node.last_in_edge = edge
                                new_reachability_probability = float(rate/total_rate) * reachability_probability
                                new_depth = depth + 1
                                
                                neighbors.append([new_node, new_reachability_probability, new_depth])

                neighbors = sorted(neighbors, key=lambda n: n[1], reverse = True)
                flag = False
                for neighbor in neighbors:
                    if neighbor[0].var_values not in visited:
                        stack.append(neighbor)
                        flag = True
                if flag:
                    current_depth = current_depth + 1
        
        for i in range(0, current_depth):
            solver.pop()
        bound = bound + 1


###################################################################
###### Greedy DFS guided by BMC and a probability threshold #######

def Guided_Greedy_DFS_probability(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound, init_threshold):
    start_time = time.time()
    bound = init_bound
    old_bound = bound
    curr_max_bound = bound
    prob_thresh = log(init_threshold)
    bp_flag = False


    nodes = {}
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)

    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()

    while True:
        current_depth = 0
        print('+++++' * 10)
        print('bound = ' + str(bound))
        print('probability threshold = ' + str(prob_thresh) + '\n')
        
        if (bound == init_bound):
            solver = Solver()
            solver.add(get_initial_state(model))

            for j in range(1, bound+1):
                solver.add(get_encoding(model, j))

            solver.push()
            x = Int(target_var + '.' + str(bound))
            property_constraint = (x==target_value)
            solver.add(property_constraint) 
            old_bound = bound 
        
        if bound == (old_bound+1):
            solver.pop()
            solver.add(get_encoding(model, bound))
            solver.push()
            x = Int(target_var + '.' + str(bound))
            property_constraint = (x==target_value)
            solver.add(property_constraint)
            old_bound = bound
        
        # elif (bound != init_bound):
        #     raise Exception("bound did not increase by one after an iteration")

        visited = {}
        stack = [[init_node, 1.0, 0]]
        last_visited_depth = 0
        
        while stack:
            node, reachability_probability, depth = stack.pop()

            if depth < current_depth:
                for i in range(0, current_depth - depth):
                    solver.pop()
                current_depth = depth
                last_visited_depth = depth
            
            if depth != 0:
                if depth != last_visited_depth:
                    solver.push()
                    x = Int('selected_reaction.' + str(depth - 1))
                    solver.add(x == node.last_in_edge.reaction)
                    last_visited_depth = depth
                    
                else:
                    solver.pop()
                    solver.push()
                    x = Int('selected_reaction.' + str(depth - 1))
                    solver.add(x == node.last_in_edge.reaction)
            
            if depth > current_depth:
                raise Exception("current_depth is behind the depth of the node at the top of the stack")

            if node.var_values not in visited:
                visited[node.var_values] = node
                if is_target(node.var_values, target_index, target_value):
                    neighbors = []
                    backtrack(diag, node, model)
                    if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                            result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                            print('prob_thresh=' + str(prob_thresh))
                            print('curr_max_bound= ' + str(curr_max_bound))
                            print('current witness information: ')
                            print('probability= ' + str(result))
                            diag_size = len(diag.edges) + len(diag.nodes)
                            print('size= ' + str(diag_size))
                            print('time= ' + str(time.time()-start_time))
                            print('='*50)
                            if time.time()-start_time > 3600:
                                quit()
                else:
                    if (diag.get_node(node.var_values)[0]):
                        backtrack(diag, node, model)
                        if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                            result = diag.check_probability(model, model_name, prism_bin, csl_prop)
                            print('prob_thresh')
                            print('curr_max_bound= ' + str(curr_max_bound))
                            print('current witness information: ')
                            print('probability= ' + str(result))
                            diag_size = len(diag.edges) + len(diag.nodes)
                            print('size= ' + str(diag_size))
                            print('time= ' + str(time.time()-start_time))
                            print('='*50)
                            if time.time()-start_time > 3600:
                                quit()
                    
                    neighbors = []
                    if (reachability_probability >= prob_thresh) and (solver.check()==sat):
                        for i, r in enumerate(model.get_reactions_vector()):
                            rate = get_reaction_rate(node.var_values, model, i)
                            total_rate = get_total_outgoing_rate(node.var_values, model)
                            if rate > 0:
                                dst_var_values = [None] * len(model.get_species_tuple())
                                for j, coefficient in enumerate(r):
                                    dst_var_values[j] = node.var_values[j] + coefficient
                                dst_var_values = tuple(dst_var_values)

                                if dst_var_values in nodes:
                                    new_node = nodes[dst_var_values]
                                else:
                                    new_node = Node()
                                    new_node.var_values = dst_var_values
                                    nodes[new_node.var_values] = new_node
                                    
                                
                                
                                edge = Edge()
                                edge.src = node
                                edge.dst = new_node
                                edge.rate = rate
                                edge.reaction = i
                                
                                new_node.add_in_edge(edge)
                                new_node.last_in_edge = edge
                                new_reachability_probability = log(float(rate/total_rate)) + reachability_probability
                                new_depth = depth + 1
                                
                                neighbors.append([new_node, new_reachability_probability, new_depth])

                neighbors = sorted(neighbors, key=lambda n: n[1], reverse = True)
                flag = False
                for neighbor in neighbors:
                    if neighbor[0].var_values not in visited:
                        stack.append(neighbor)
                        flag = True
                if flag:
                    current_depth = current_depth + 1
        
        for i in range(0, current_depth):
            solver.pop()
        if bp_flag:
            bound = bound + 1
            if bound > curr_max_bound:
                curr_max_bound = bound
                bp_flag = False
        else:
            prob_thresh = prob_thresh - 1
            bound = init_bound
            bp_flag = True