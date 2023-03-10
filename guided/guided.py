import numpy as np
import Models
from z3 import *
from Graph import graph, node, edge
from Utils import *


import os, sys, json, time

f = open(sys.argv[1])
json_data = json.load(f)
model_name = json_data['model']
starting_bound = int(json_data['starting_bound'])
csl_prop = json_data['csl_property']
prism = json_data['prism_binary']
prob_bound = float(json_data['probability_bound'])
property_var = json_data['property_variable']
property_val = int(json_data['property_value'])
A5 = json_data['A5']
reaction_subset = list(A5.split(","))

#input: a state_vector(variable assignments) and a bound - ouput: True/False (is target achievable from
#this state within 'bound' number of transitions?)
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

#input: a node - ouput: True/False   
def is_target(node, property_var, property_val):
    if node.var_dict[property_var] == property_val:
        return True
    return False

#input a node and a graph object. Recursive function which follows the parent
#link in node and expands the trace until it reaches the initial state.
def construct_trace(node):
    trace = graph()
    while(node.parent != None):
        reaction = node.taken_reaction
        n1 = node.parent
        n2 = node
        edge_ = edge(n1, n2, reaction)
        trace.add_edge(edge_)
        node = node.parent
    return trace

#input: var_dict of a node - output: a dictionary contatining transitions with p>0
def return_r_dict(model, var_dict):
    r_dict = {}
    #print (var_dict)
    total_rate = 0
    for i, r in enumerate(model.reactions_dict()):
        comb = 1
        rate = model.reaction_rates()[i]
        lhs = model.reactions_dict()[r]
        for j, v in enumerate(lhs[0]):
            s = model.species_vector()[j]
            for c in range(v): 
                comb = comb * var_dict[s]
        rate = rate * comb
        total_rate = total_rate + rate

    r_rate = 0
	# 2R1 + R2 --K--> R3. the rate would be: #(R1) * #(R1) * #(R2) * K
    for i, r in enumerate(model.reactions_dict()):
        comb = 1
        rate = model.reaction_rates()[i]
        lhs = model.reactions_dict()[r]
        for j, v in enumerate(lhs[0]):
            s = model.species_vector()[j]
            for c in range(v): 
                comb = comb * var_dict[s]
        rate = rate * comb
        rate = float(float(rate) / float(total_rate))
        if rate!=0:
            r_dict[r] = rate 
    return r_dict

def find_trace(model, curr_node, bound,  property_val, property_var, probability = 0):
    # print('++++')
    # print("current node var dict:")
    # print(curr_node.var_dict)
    # print("------")
    # print("current node outgoing transitions:")
    # print(curr_node.r_dict)
    # print("------")
    # print('++++')
    #case1: target is reached
    if is_target(curr_node, property_val=property_val, property_var=property_var):
        #trace would be a graph object contatining target node
        trace = graph()
        trace = construct_trace(curr_node)
        taken_reaction = curr_node.taken_reaction
        prev_node = curr_node.parent
        del prev_node.r_dict[taken_reaction]
        return trace, prev_node, bound+1, False
        
        
        
            
    #case2: outgoing transitions from this node are exhausted
    elif len(curr_node.r_dict)==0:
        parent_node = curr_node.parent
        if parent_node==None:
            print("reached init state")
            return None, None, None, True
        parent_bound = bound + 1
        curr_node.child = None
        del parent_node.r_dict[curr_node.taken_reaction]
        return find_trace(model=model, curr_node=parent_node, bound=parent_bound, property_val=property_val, property_var=property_var)

    #case3: there are some outgoing transitions from this node
    else:
        max_prob = 0
        max_prob_key = None
        
        for e in curr_node.r_dict:
            if curr_node.r_dict[e] > max_prob:
                max_prob = curr_node.r_dict[e]
                max_prob_key = e
        curr_state_vector = [None] * len(model.species_vector())
        
        for i, e in enumerate(curr_node.var_dict):
            curr_state_vector[i] = curr_node.var_dict[e]
        next_state_vector = (np.array(curr_state_vector) + np.array(model.reactions_dict()[max_prob_key][2])).tolist()
        
        #case3a: if selecting the highest prob transition wouldn't work delete that transition and check again
        if not(check_sat(model, next_state_vector, bound-1, target_value=property_val, target_var=property_var)):
            del curr_node.r_dict[max_prob_key]
            if not(curr_node.child == None):
                if curr_node.child.taken_reaction == max_prob_key:
                    curr_node.child = None
            return find_trace(model, curr_node, bound, property_val=property_val, property_var=property_var)
        
        #case3b: move the head to the new state
        else:        
            var_dict_next ={}
            for i, s in enumerate(model.species_vector()):
                var_dict_next[s] = next_state_vector[i]
            
            if (curr_node.child == None):
                next_node = node(var_dict=var_dict_next, parent=curr_node, r_dict=return_r_dict(model, var_dict_next), taken_reaction = max_prob_key)
                curr_node.child = next_node
            elif (curr_node.child.var_dict == var_dict_next):
                next_node = curr_node.child
            else:
                next_node = node(var_dict=var_dict_next, parent=curr_node, r_dict=return_r_dict(model, var_dict_next), taken_reaction = max_prob_key)
                curr_node.child = next_node
        
            return find_trace(model, next_node, bound-1, property_val=property_val, property_var=property_var)

########################################################################
constructor = getattr(Models, model_name)
model = constructor()

#finding the length of the shortest path to target state
min_bound = 0
while (True):
    if check_sat(model, model.initial_state(), min_bound, target_var=property_var, target_value=property_val):
        break
    else:
        min_bound = min_bound + 1

#print (min_bound)
#pointing the current node to the initial state of the graph
initial_state_vector = model.initial_state()
var_dict = {}
for i, s in enumerate(model.species_vector()):
    var_dict[s] = initial_state_vector[i]


returned_graph = graph()

min_bound_copy = min_bound

terminate_flag = False
with open('./results/' + model_name + '/' + model_name + '.results', mode = 'w', encoding= 'ascii') as f:
    f.truncate()
    start_time = time.time()
    j = 0
    while(not(terminate_flag)):
        print(j)
        if j<6:

            node_ = node(var_dict, r_dict = return_r_dict(model, var_dict))
            node_.make_initial()
            curr_node = node_
            min_bound = min_bound_copy + j
            j = j+1
            
            for i in range(100):
                a, b, bound_, terminate = find_trace(model, curr_node, min_bound, property_var=property_var, property_val=property_val)
                if terminate:
                    print("exhausted")
                    break

                for i, e in enumerate(a.edge_list):
                    returned_graph.add_edge(e)
                    #print(str(i+1) + ":" + str(e.n1.var_dict) + " ->" + str(e.n2.var_dict) + "|" + str(e.reaction))
                # prob = returned_graph.model_check(model, model_name, prism, csl_prop)
                # print(prob)
                # print ("=======")
                curr_node = b
                min_bound = bound_
                a = None

        scaffold(returned_graph, model, 20, 500, property_var, property_val, reaction_subset)
        prob = returned_graph.model_check(model, model_name, prism, csl_prop)
        elapsed_time = time.time()
        if (elapsed_time-start_time)>1800:
            terminate_flag = True
        print(prob)
        print(elapsed_time)
        f.write('# of nodes: ' + str(len(returned_graph.node_list)))
        f.write('\n')
        f.write('# of edges: ' + str(len(returned_graph.edge_list)))
        f.write('\n')
        f.write('probability = ' + str(prob))
        f.write('\n')
        f.write('elapsed time: ' + str(elapsed_time-start_time))
        f.write('\n')
        f.write('='*40)
        f.write('\n')
    f.close()
















####### until here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#print("seed path found")

#Yeast
# reaction_subset = [1,1,1,1,1,1,1,1]
# while (True):
#     scaffold(trace, model, 20, 200, "Gbg", 50, reaction_subset)
#     print(trace.model_check(model, "yeast_polarization", prism, "P=? [true U<=20 Gbg=50]"))





###################################################
# test graph

###############
# n0 = node(var_dict={"S1":1, "S2":2, "S3":2})
# n1 = node(var_dict={"S1":1, "S2":3, "S3":3})
# n2 = node(var_dict={"S1":1, "S2":1, "S3":1})
# n2.make_initial()
# n3 = node(var_dict={"S1":1, "S2":2, "S3":3})
# n4 = node(var_dict={"S1":1, "S2":1, "S3":2})
# n5 = node(var_dict={"S1":1, "S2":1, "S3":3})
# n6 = node(var_dict={"S1":-1, "S2":-1, "S3":-1})

# e1=edge(n1=n2, n2=n0, reaction=1)
# e2=edge(n1=n2, n2=n4, reaction=2)
# e3=edge(n1=n0, n2=n1, reaction=1)
# e4=edge(n1=n0, n2=n3, reaction=2)
# e5=edge(n1=n4, n2=n3, reaction=1)
# e6=edge(n1=n4, n2=n5, reaction=2)
# e7=edge(n1=n1, n2=n6, reaction=1)
# e8=edge(n1=n1, n2=n6, reaction=2)
# e9=edge(n1=n3, n2=n6, reaction=1)
# e10=edge(n1=n3, n2=n6, reaction=2)
# e11=edge(n1=n5, n2=n6, reaction=1)
# e12=edge(n1=n5, n2=n6, reaction=2)

# test_graph = graph()
# test_graph.add_edge(e1)
# test_graph.add_edge(e2)
# test_graph.add_edge(e3)
# test_graph.add_edge(e4)
# test_graph.add_edge(e5)
# test_graph.add_edge(e6)
# test_graph.add_edge(e7)
# test_graph.add_edge(e8)
# test_graph.add_edge(e9)
# test_graph.add_edge(e10)
# test_graph.add_edge(e11)
# test_graph.add_edge(e12)

# prob = test_graph.model_check(model, model_name, prism, csl_prop)
# print(prob)



# print("="*60)
###############
# test graph
###################################################
