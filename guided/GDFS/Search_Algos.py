from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, check_sat
from Graph_ import Node, Edge, Graph
def construct_trace(graph, node):
    curr_node = node
    flag = True
    while((not graph.get_node(curr_node.var_values)[0]) or flag):
        flag = False
        graph.add_edge(curr_node.incoming_edge)
        curr_node = curr_node.incoming_edge.src

#guided greedy DFS where all the nodes allowed by the minmax dictionary search stack. 
#All nodes added to the search stack are added to the graph.
def ggdfs(graph, start_node, model, target_index, target_value, min_max):
    flag = False
    visited = {}
    stack = [start_node]
    while stack:
        node = stack.pop()
        #print(node.get_var_values())
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

#guided greedy DFS where all the nodes allowed by the minmax dictionary search stack. 
#Whenever the search reaches a target or a node already in the graph,
#the algorithm traces back and adds all the edges to the graph. the size of the returned graph would
#be smaller than ggdfs
def ggdfs_trace_back(graph, start_node, model, target_index, target_value, min_max):
    flag = False
    visited = {}
    stack = [start_node]
    while stack:
        node = stack.pop()
        if (graph.get_node(node.var_values)[0] or is_target(node.var_values, target_index, target_value)):
            if (not node.initial_state):
                construct_trace(graph, node)
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
                        for j, coefficient in enumerate(r[2]):
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
                            new_node.incoming_edge = edge
                            neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append(neighbor)
    return flag

#guided greedy DFS where all the nodes with reachability probability greater than the threshold 
#are added to the search stack. All nodes added to the search stack are added to the graph.
def ggdfs_prob(graph, start_node, model, target_index, target_value, prob_thresh, min_max):
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
                        for j, coefficient in enumerate(r[2]):
                            dst_var_values[j] = node.var_values[j] + coefficient
                        dst_var_values = tuple(dst_var_values)
                        if check_sat(dst_var_values, min_max):
                            new_node = Node()
                            new_node.var_values = dst_var_values
                            if graph.get_node(dst_var_values)[0]:
                                new_node = graph.get_node(dst_var_values)[1]
                            new_node.reachability_probability = new_node.reachability_probability + float(rate/total_rate)*node.reachability_probability

                            edge = Edge()
                            edge.src = node
                            edge.dst = new_node
                            edge.rate = rate
                            edge.reaction = i

                            if not graph.get_edge(edge.get_tuple())[0]:
                                graph.add_edge(edge)
                                # node.reachability_probability = node.reachability_probability + float(rate/total_rate)
                
                            if new_node.reachability_probability>= prob_thresh:
                                neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append(neighbor)
    return flag

#guided Greedy DFS where all the nodes with reachability probability greater than the threshold 
#are added to the search stack. Whenever the search reaches a target or a node already in the graph,
#the algorithm traces back and adds all the edges to the graph. the size of the returned graph would
#be smaller than ggdfs_prob
def ggdfs_prob_trace_back(graph, start_node, model, target_index, target_value, prob_thresh, min_max):
    flag = False
    visited = {}
    stack = [start_node]
    while stack:
        node = stack.pop()
        if (graph.get_node(node.var_values)[0] or is_target(node.var_values, target_index, target_value)):
            if (not node.initial_state):
                construct_trace(graph, node)
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
                        for j, coefficient in enumerate(r[2]):
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

                            if new_node.reachability_probability> prob_thresh:
                                new_node.incoming_edge = edge
                                neighbors.append(new_node) 
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append(neighbor)
    return flag