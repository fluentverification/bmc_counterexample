from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, check_sat
from Graph import Node, Edge, Graph
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
def ggdfs_prob(graph, start_node, model, target_index, target_value, min_max, prob_thresh, size):
    flag = False
    visited = {}
    stack = [(start_node, 1.0)]
    while stack:
        node, prob = stack.pop()
        
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
                        if check_sat(dst_var_values, min_max) and prob>prob_thresh:
                            new_node = Node()
                            new_node.var_values = dst_var_values
                            new_node.reachability_probability = (float(rate/total_rate)*prob)
                            edge = Edge()
                            edge.src = node
                            edge.dst = new_node
                            edge.rate = rate
                            edge.reaction = i

                            if prob>prob_thresh:
                                if new_node.var_values in visited:
                                    old_prob = visited[new_node.var_values].reachability_probability
                                    if new_node.reachability_probability > old_prob:
                                        del visited[new_node.var_values]

                                # if new_node.var_values in visited:
                                #     del visited[new_node.var_values]
                                if not graph.get_edge(edge.get_tuple())[0]:
                                    graph.add_edge(edge)
                                neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append((neighbor, new_node.reachability_probability))
    return flag



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
                        if check_sat(dst_var_values, min_max) and ((prob*float(rate/total_rate))>prob_thresh):
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
