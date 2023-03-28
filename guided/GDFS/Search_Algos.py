from Utils import is_target, get_reaction_rate, get_total_outgoing_rate, check_sat
from Graph_ import Node, Edge, Graph

def gdfs_prob(graph, start_node, model, prism_bin, csl_prop, target_index, target_value, prob_thresh, start_time, min_max):
    flag = False
    graph_size = len(graph.nodes) + len(graph.edges)
    
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
                        for j, coefficient in enumerate(r[2]):
                            dst_var_values[j] = node.var_values[j] + coefficient
                        dst_var_values = tuple(dst_var_values)
                        if check_sat(dst_var_values, min_max):
                            bound = node.bound
                            new_node = Node()
                            new_node.var_values = dst_var_values
                            new_node.bound = bound-1
                            
                            new_node.reachability_probability = float(rate/total_rate)*node.reachability_probability
                            
                            edge = Edge()
                            edge.src = node
                            edge.dst = new_node
                            edge.rate = rate
                            edge.reaction = i

                            if not graph.get_edge(edge.get_tuple())[0]:
                                graph.add_edge(edge)
                                node.reachability_probability = node.reachability_probability + float(rate/total_rate)
                
                            if new_node.reachability_probability> prob_thresh:
                                neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.reachability_probability, reverse = True)
            for neighbor in neighbors:
                if neighbor.var_values not in visited:
                    stack.append(neighbor)
            #print(len(graph.nodes) + len(graph.edges))
            if len(graph.nodes) + len(graph.edges) - graph_size > 20000:
                graph_size = len(graph.nodes) + len(graph.edges)
                current_time = time.time()
                with open('result' + '.txt', mode='a', encoding='ascii') as f:
                    prob = graph.check_probability(model, 'test', prism_bin, csl_prop)
                    size = len(graph.nodes) + len(graph.edges)
                    f.write('prob: ' + str(prob))
                    f.write("\n")
                    f.write('size: ' + str(size))
                    f.write("\n")
                    f.write('time: ' + str(current_time-start_time))
                    f.write("\n")
                    f.write("="*50)
                    f.write("\n")
                    print(prob)
                    print(size)
                    print('='*50)
                    f.close()

                if (current_time - start_time > 440000):
                    break

    return flag
