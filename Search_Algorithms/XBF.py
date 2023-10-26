from heapq import heapify, heappush, heappop
from math import log
import time
from Parser import Parser
from Graph import Node, Edge, Graph, check_probability
from Utils import get_total_outgoing_rate, get_reaction_rate, is_target

def CEX_GEN(json_data):
    start_time = time.time()

    #parsing parameters from json elements
    model_path = json_data['model_path']
    model_name = json_data['model_name']
    # K = int(json_data['starting_bound'])
    csl_prop_lb = json_data['csl_property']
    # csl_prop_ub = json_data['csl_property_upper_bound']
    prism_bin = json_data['prism_binary']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    mc_step = int(json_data['model_check_step'])

    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    print('Starting [XBF]... \n')
    XBF(model, target_index, target_value, mc_step, model_name, prism_bin, csl_prop_lb, start_time)


def h(reach_prob):
    return 0-reach_prob

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

def XBF(model, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, start_time):
    init_node = Node()
    init_node.initial_state = True
    init_node.var_values = model.get_initial_state()
    init_node.reachability_probability = 0.0
    nodes = {}
    nodes[init_node.var_values] = init_node
    # diag <-- empty graph
    # closed <-- empty hash table
    # open_ <-- empty priority queue (implemented as a heap)
    # map_ <-- mapping node's var_values to a list of size 2 containing the hueristic value and the node object
    # map_ gives access to the list elements of open_. after an element of the list of size 2 is changed, the 
    # open_ should be heapified again in order to maintain the logic of the heap
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)
    open_ = []
    closed = {}
    map_ = {}

    node_list = [h(0.0), init_node]
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
            # dst_var_values = [None] * len(model.get_species_tuple())
            # for j, coefficient in enumerate(r):
            #     dst_var_values[j] = s.var_values[j] + coefficient
            #

            if rate > 0 : #and check_sat(dst_var_values, min_max_list):
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
                if (s.reachability_probability + log((rate/total_rate), 10)) > s_prime.reachability_probability:
                    s_prime.reachability_probability = s.reachability_probability + log((rate/total_rate), 10)

                new_h = h(s_prime.reachability_probability)

                if (is_target(s_prime.var_values, target_index, target_value)
                    or diag.get_node(s_prime.var_values)[0]):
                    backtrack(diag, s_prime, model)
                    if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
                        result = check_probability(diag, model, model_name, prism_bin, csl_prop)
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
                            heapify(open_)
                    elif (s_prime.var_values in closed):
                        if new_h < closed[s_prime.var_values][0]:
                            del closed[s_prime.var_values]
                            node_list = [new_h, s_prime]
                            heappush(open_, node_list)
                            if s_prime.var_values in map_:
                                raise Exception("a node is both in open and closed")
                            map_[s_prime.var_values] = node_list