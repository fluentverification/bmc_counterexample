############################################################
############################################################
#Notes: 
#	1-this module works only for models with a single initial state
#	2-only properties of type "a U<=t b" where b is of the form "var=value"
#	  is accepted currently	
############################################################
############################################################

from z3 import *
import Models
from Graph import graph, node
from Utils import *
import sys
import json
import numpy as np
import os
import time
############################################################
############################################################
############################################################

#read json elements
f = open(sys.argv[1])
json_data = json.load(f)
model_name = json_data['model']
starting_bound = int(json_data['starting_bound'])
csl_prop = json_data['csl_property']
csl_prop_sink = json_data['csl_property_sink']
prism = json_data['prism_binary']
prob_bound = float(json_data['probability_bound'])
property_var = json_data['property_variable']
property_val = int(json_data['property_value'])
mc_step = int(json_data['model_check_step'])
scaffold_bound_limit = int(json_data['scaffold_bound_limit'])
scaffold_count_limit = int(json_data['scaffold_count_limit'])
start_scaffold = int(json_data['start_scaffold'])
model_bound_prob_threshold = float(json_data['model_bound_prob_threshold'])
A1 = json_data['A1']
A2 = json_data['A2']
A3 = json_data['A3']
A4 = json_data['A4']
A5 = json_data['A5']
A1 = list(A1.split(","))
A2 = list(A2.split(","))
A3 = list(A3.split(","))
A4 = list(A4.split(","))
A5 = list(A5.split(","))

#counterexample will be saved in results folder
if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)

#'model' will point to the specific model class within 
#models module
constructor = getattr(Models, model_name)
model = constructor()
graph = graph()

#adding the initial state to the graph
initial_state_vector = model.initial_state()
var_dict = {}
for i, s in enumerate(model.species_vector()):
    var_dict[s] = initial_state_vector[i]
node_ = node(var_dict)
node_.make_initial()
graph.add_node(node_)

#
count = 0
prob = 0
prob_prev = 0
prob_sink = 0
model_bound = 0
flag = True
start_time = time.time()
excluded_paths = []

with open('./results/' + model_name + '/' + model_name + '.results', mode = 'w', encoding= 'ascii') as f:
    f.truncate()
	
    L = 19
    U = 19
    
    solver = Solver()
    solver.add(get_initial_state(model))
    solver.add(Real("rate.0")==1.0)
    solver.add(get_encoding_LU(model, L, U, property_var, property_val))


    while prob <= prob_bound:

        while (solver.check() == sat):
            count = count + 1
            path = solver.model()
            print(path[Real("rate.18")])
            ep = exclude_path(path)
            solver.add(ep)
            graph.add_path(path)

            # if (count % start_scaffold == 0):
            #     scaffold(graph, model, bound-1, scaffold_count_limit, property_var, property_val, reaction_subset)

            prob = graph.model_check(model, model_name, prism, csl_prop)
            prob_sink = graph.model_check(model, model_name, prism, csl_prop_sink)
            
            elapsed_time = time.time()
            print('# of nodes: ' + str(len(graph.node_list)))
            print('# of edges: ' + str(len(graph.edge_list)))
            print('probability = ' + str(prob))
            print('upper-bound probability = ' + str(prob_sink))
            print('elapsed time: ' + str(elapsed_time-start_time))
            print('='*40)
            f.write('# of nodes: ' + str(len(graph.node_list)))
            f.write('\n')
            f.write('# of edges: ' + str(len(graph.edge_list)))
            f.write('\n')
            f.write('probability = ' + str(prob))
            f.write('\n')
            f.write('upper-bound probability = ' + str(prob_sink+ prob))
            f.write('\n')
            f.write('elapsed time: ' + str(elapsed_time-start_time))
            f.write('\n')
            f.write('='*40)
            f.write('\n')


        elapsed_time = time.time()
        if (elapsed_time-start_time>1800):
            break
        
    f.close()

f.close()
