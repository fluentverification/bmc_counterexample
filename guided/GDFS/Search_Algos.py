from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, check_sat
from Graph import Node, Edge, Graph
from queue import PriorityQueue
from heapq import heapify, heappush, heappop

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

def h(reach_prob, var_value, model):
    return reach_prob

def backtrack(graph, node, model):
    pass

def XBF(model, target_index, target_value, mc_step):
    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()
    init_node.reachability_probability = 1.0 
    # diag <-- empty graph
    # closed <-- empty hash table
    # open_ <-- empty priority queue (implemented using heap functions)
    # map_ <-- mapping node's var_values to a list of size 2 containing the hueristic value and the node object
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)
    open_ = []
    closed = {}
    map_ = {}

    node_list = [h(1.0, init_node.var_values, model), init_node]
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
            if rate > 0:
                dst_var_values = [None] * len(model.get_species_tuple())
                for j, coefficient in enumerate(r):
                    dst_var_values[j] = s.var_values[j] + coefficient
                dst_var_values = tuple(dst_var_values)
                s_prime = Node()
                s_prime.var_values = dst_var_values
                edge = Edge()
                edge.src = s
                edge.dst = s_prime
                edge.rate = rate
                edge.reaction = i
                s_prime.add_in_edge(edge)
                s_prime.reachability_probability = s.reachability_probability * float(rate/total_rate)

                if is_target(s_prime.var_values, target_index, target_value):
                    backtrack(diag, s_prime, model)
            
                new_h = h(s_prime.reachability_probability, s_prime.var_values, model)

                if (s_prime.var_values not in map_ 
                    and s_prime.var_values not in closed):
                    node_list = [new_h, s_prime]
                    heappush(open_, node_list)
                    map_[s_prime.var_values] = node_list
                else:
                    if (s_prime.var_values in map_):
                        if new_h < map_[s_prime.var_values][0]:
                            map_[s_prime.var_values][0] = new_h
                            heapify(open_)
                    elif (s_prime.var_values in closed):
                        if new_h < closed[s_prime.var_values][0]:
                            del closed[s_prime.var_values]
                            node_list = [new_h, s_prime]
                            heappush(open_, node_list)
                            map_[s_prime.var_values] = node_list


