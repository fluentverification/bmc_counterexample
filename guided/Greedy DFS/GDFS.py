from Utils import get_reaction_rate, get_total_outgoing_rate, set_
from Graph_ import Node, Edge, Graph
from Models import enzymatic_futile_cycle, test, motility_regulation, yeast_polarization
from z3 import *
import time
import sys

def check_sat(model, state_vector, bound, target_index, target_value):
    vars = []
    for i,r in enumerate (model.reactions_vector()):
        x = Int("n" + str(i))
        vars.append(x)
    
    constraints = []
    #first constraint (n1,n2,...>=0)
    for i in vars:
        constraints.append(i>=0)

    #second constraint (n1+n2+...=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum<=bound))
    
    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.reactions_vector()):
        if r[2][target_index]!=0:
            vars.append([Int("n" + str(i)), r[2][target_index]])
    sum = state_vector[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    #fourth constraint (species population >=0)
    for i, s in enumerate(model.species_vector()):
        if i==target_index:
            continue
        vars = []
        for j, r in enumerate(model.reactions_vector()):
            if r[2][i]!=0:
                vars.append([Int("n" + str(j)), r[2][i]])
        sum = state_vector[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0) 

    solver = Solver()
    solver.add(And(constraints))
    if (solver.check()==sat):
        return True
    else:
        return False

def is_target(var_values, target_index, target_value):
    if var_values[target_index] == target_value:
        return True
    return False


def gdfs(graph, start_node, model, prism_bin, csl_prop, target_index, target_value):
    start_time = time.time()
    graph_size = 0
    
    visited = set_()
    stack = [start_node]
    while stack:
        node = stack.pop()
        #print(node.get_var_values())
        if not node in visited:
            visited.add(node)
            if is_target(node.get_var_values(), target_index, target_value):
                neighbors = []
            else:
                neighbors = []
                for i, r in enumerate(model.reactions_vector()):
                    rate = get_reaction_rate(node.get_var_values(), model, i)
                    total_rate = get_total_outgoing_rate(node.get_var_values(), model)
                    if rate > 0:
                        dst_var_values = [None] * len(model.species_vector())
                        for j, coefficient in enumerate(r[2]):
                            dst_var_values[j] = node.get_var_values()[j] + coefficient
                        if check_sat(model, dst_var_values, node.get_bound(), target_index, target_value):
                            bound = node.get_bound()
                            new_node = Node(var_values = dst_var_values, bound = bound-1)
                            if not new_node in visited:
                                new_node.set_reachability_probability(float(rate/total_rate))
                                graph.add_node(new_node)
                            else:
                                for n in visited:
                                    if n.equals(new_node):
                                        new_node = n
                                        if new_node.get_bound() < (bound - 1):
                                            new_node.set_bound(bound - 1)
                            edge = Edge()
                            edge.set_src(node)
                            edge.set_dst(new_node)
                            edge.set_rate(rate)
                            edge.set_prob_dtmc(float(rate/total_rate))
                            node.add_outgoing_edge(edge)
                            graph.add_edge(edge)
                            neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.get_reachability_probability(), reverse = True)
            for neighbor in neighbors:
                if neighbor not in visited:
                    stack.append(neighbor)
            if graph.get_size() - graph_size > 2000:
                graph_size = graph.get_size()
                current_time = time.time()
                with open(model.name() + '.txt', mode='a', encoding='ascii') as f:
                    prob = graph.model_check(model, prism_bin, csl_prop)
                    size = graph.get_size()
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
        current_time = time.time()
        if (current_time - start_time > 3600):
            break

                    
def gdfs_prob(graph, start_node, model, prism_bin, csl_prop, target_index, target_value, prob_thresh, start_time):
    graph_size = graph.get_size()
    
    visited = set_()
    stack = [start_node]
    while stack:
        node = stack.pop()
        #print(node.get_var_values())
        if not node in visited:
            visited.add(node)
            if is_target(node.get_var_values(), target_index, target_value):
                neighbors = []
            else:
                neighbors = []
                for i, r in enumerate(model.reactions_vector()):
                    rate = get_reaction_rate(node.get_var_values(), model, i)
                    total_rate = get_total_outgoing_rate(node.get_var_values(), model)
                    if rate > 0:
                        dst_var_values = [None] * len(model.species_vector())
                        for j, coefficient in enumerate(r[2]):
                            dst_var_values[j] = node.get_var_values()[j] + coefficient
                        if check_sat(model, dst_var_values, node.get_bound(), target_index, target_value):
                            bound = node.get_bound()
                            new_node = Node(var_values = dst_var_values, bound = bound-1)
                            
                            new_node.set_reachability_probability(float(rate/total_rate)*node.get_reachability_probability())
                            contain_node = graph.contain_node(new_node)
                            if not contain_node[0]:
                                graph.add_node(new_node)
                            else:
                                new_node = contain_node[1]
                                edge = Edge()
                                edge.set_src(node)
                                edge.set_dst(new_node)
                                edge.set_rate(rate)
                                edge.set_prob_dtmc(float(rate/total_rate))
                                edge.set_reaction(i)
                                if not graph.contain_edge(edge)[0]:
                                    new_node.set_reachability_probability(new_node.get_reachability_probability()+float(rate/total_rate))
                                if bound - 1 > new_node.get_bound():
                                    new_node.set_bound(bound - 1)      
                            edge = Edge()
                            edge.set_src(node)
                            edge.set_dst(new_node)
                            edge.set_rate(rate)
                            edge.set_prob_dtmc(float(rate/total_rate))
                            edge.set_reaction(i)
                            if not graph.contain_edge(edge)[0]:
                                graph.add_edge(edge)
                            if not node.contain_outgoing_edge(edge):
                                node.add_outgoing_edge(edge)

                            if new_node.get_reachability_probability()> prob_thresh:
                                neighbors.append(new_node)
            neighbors = sorted(neighbors, key=lambda n: n.get_reachability_probability(), reverse = True)
            for neighbor in neighbors:
                if neighbor not in visited:
                    stack.append(neighbor)
            if graph.get_size() - graph_size > 4000:
                graph_size = graph.get_size()
                current_time = time.time()
                with open(model.name() + '.txt', mode='a', encoding='ascii') as f:
                    prob = graph.model_check(model, prism_bin, csl_prop)
                    size = graph.get_size()
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




prism_bin = "../../../prism-4.7-src/prism/bin/prism"
sys.setrecursionlimit(1000000)
#Enzymatic Futile cycle
# csl_prop = "P=? [true U<=100 S5=40]"
# target_index = 4
# target_value = 40
# #original bound is 19
# bound = 19
# model = enzymatic_futile_cycle()
# prob_thresh = 1E-29

#Yeast Polarization

csl_prop = "P=? [true U<=20 Gbg=50]"
target_index = 5
target_value = 50
#original bound is 100
bound = 100
model = yeast_polarization()
prob_thresh = 1E-2

#Motility Regulation
# csl_prop = "P=? [true U<=10 CodY=19]"
# target_index = 6
# target_value = 19
# #original bound is 9
# bound = 9
# model = motility_regulation()
# prob_thresh = 1E-10

start_time = time.time()
init_node = Node()
init_node.set_is_initial(True)
init_node.set_var_values(model.initial_state())
init_node.set_bound(bound)
init_node.set_reachability_probability(1.0)
graph = Graph()
graph.add_node(init_node)
#gdfs(graph=graph, start_node=init_node, model=model, prism_bin=prism_bin, csl_prop=csl_prop, target_index=target_index, target_value=target_value)




while (True):
    gdfs_prob(graph=graph, start_node=init_node, model=model, prism_bin=prism_bin, csl_prop=csl_prop, target_index=target_index, target_value=target_value, prob_thresh=prob_thresh, start_time=start_time)
    prob_thresh = prob_thresh / 2.0
    init_node.set_bound(init_node.get_bound() + 1)
    print(prob_thresh)
    current_time = time.time()
    if (current_time - start_time > 440000):
        break











# def dfs_greedy_iterative(graph, start_node):
#     visited = set()
#     stack = [start_node]
#     while stack:
#         node = stack.pop()
#         if node not in visited:
#             visited.add(node)
#             neighbors = graph[node]
#             neighbors = sorted(neighbors, key=lambda x: graph[node][x], reverse=True)
#             for neighbor in neighbors:
#                 if neighbor not in visited:
#                     stack.append(neighbor)
#     return visited