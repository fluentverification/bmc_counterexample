import numpy as np
import Models
from z3 import *
from Graph import graph, node, edge


def check_sat(model, state_vector, bound, target_var, target_value):
    vars = []
    for r in model.reactions_dict():
        x = Int("n" + str(r))
        vars.append(x)
    
    constraints = []
    #first constraint (n1,n2,...>=0)
    for i in vars:
        constraints.append(i>=0)

    #second constraint (n1+n2+...=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum==bound))
    
    #third constraint (reaching the target)
    vars = []
    for i, s in enumerate(model.species_vector()):
        if s==target_var:
            target_index = i
    for i, r in enumerate(model.reactions_dict()):
        if model.reactions_dict()[r][2][target_index]!=0:
            vars.append([Int("n" + str(r)), model.reactions_dict()[r][2][target_index]])
    sum = state_vector[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    #fourth constraint (species population >=0)
    for i, s in enumerate(model.species_vector()):
        if s==target_var:
            continue
        vars = []
        for j, r in enumerate(model.reactions_dict()):
            if model.reactions_dict()[r][2][i]!=0:
                vars.append([Int("n" + str(r)), model.reactions_dict()[r][2][i]])
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
    
def is_target(node):
    if node.var_dict["Gbg"] == 50:
        return True
    return False

def construct_trace(node, trace):
    print(node.var_dict)
    if node.parent == None:
        return
    
    n1 = node.parent
    n2 = node
    edge_ = edge(n1, n2)
    trace.add_edge(edge_)
    construct_trace(n1, trace)
########################################################################


constructor = getattr(Models, "yeast_polarization")
model = constructor()

#finding the length of the shortest path to target state
min_bound = 0
while (True):
    if check_sat(model, model.initial_state(), min_bound, target_var="Gbg", target_value=50):
        break
    else:
        min_bound = min_bound + 1

#pointing the current node to the initial state of the graph
initial_state_vector = model.initial_state()
var_dict = {}
for i, s in enumerate(model.species_vector()):
    var_dict[s] = initial_state_vector[i]
node_ = node(var_dict)
node_.make_initial()
curr_node = node_


#select the next node by selecting the highest (possible) transition
#if you reach target, trace back parent nodes to construct the trace
while (True):
    if is_target(curr_node):
        #trace would be a graph object contatining target node
        trace = graph()
        construct_trace(curr_node, trace)
        break

#

# var_dict_2 = {"R": 2, "L": 2, "RL": 2, "G": 2, "Ga": 2, "Gbg": 2,"Gd": 2}
# node_2 = node(var_dict_2, parent = node_)

# var_dict_3 = {"R": 3, "L": 3, "RL": 3, "G": 3, "Ga": 3, "Gbg": 3,"Gd": 3}
# node_3 = node(var_dict_3, parent = node_2)

# var_dict_4 = {"R": 4, "L": 4, "RL": 4, "G": 4, "Ga": 4, "Gbg": 4,"Gd": 4}
# node_4 = node(var_dict_4, parent = node_3)

# trace = graph()
# construct_trace(node_4, trace)
# trace.to_file("test", model, "yeast_polarization")