############################################################
############################################################
#Notes: 
#	1-this module works only for models with a single initial state
#	2-only properties of type "a U<=t b" where b is of the form "var=value"
#	  is accepted currently	
############################################################
############################################################

from z3 import *
import models
from Graph import graph, node
from utils import *
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
prism = json_data['prism_binary']
prob_bound = float(json_data['probability_bound'])
property_var = json_data['property_variable']
property_val = int(json_data['property_value'])
mc_step = int(json_data['model_check_step'])
cp_bound = int(json_data['construct_path_bound'])
base = int(json_data['base'])
step = int(json_data['step'])
#counterexample will be saved in results folder
if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)

#'model' will point to the specific model class within 
#models module
constructor = getattr(models, model_name)
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
start_time = time.time()

with open('./results/' + model_name + '/' + model_name + '.results', mode = 'w', encoding= 'ascii') as f:
    f.truncate()
	
    bound = starting_bound
    while prob <= prob_bound:
		 
        solver = Solver()
        solver.add(get_initial_state(model))

        for j in range(1, bound+1):
            solver.add(get_encoding(model, j))

        x = Int(property_var + '.' + str(bound))
        property_constraint = (x==property_val)
			
        while (solver.check(property_constraint) == sat):
            path = solver.model()
            solver.add(exclude_path(path))
            graph.add_path_o(path)
            prob = graph.model_check(model, model_name, prism, csl_prop)
            print('# of nodes: ' + str(len(graph.node_list)))
            print('# of edges: ' + str(len(graph.edge_list)))
            print('probability = ' + str(prob))
            print('='*40)
            count = count + 1
            # if count>2:
            #     count = 0
            #     count, prob1, terminate = construct_path_o(graph, model, model_name, prism, csl_prop, cp_bound, count, prob, prob_bound, property_var, property_val, mc_step)
            #     prob = graph.model_check(model, model_name, prism, csl_prop)
            #     elapsed_time = time.time()
            #     # print('# of nodes: ' + str(len(graph.node_list)))
            #     # print('# of edges: ' + str(len(graph.edge_list)))
            #     # print('probability = ' + str(prob))
            #     # print('elapsed time: ' + str(elapsed_time-start_time))
            #     # print('='*40)
            #     f.write('# of nodes: ' + str(len(graph.node_list)))
            #     f.write('\n')
            #     f.write('# of edges: ' + str(len(graph.edge_list)))
            #     f.write('\n')
            #     f.write('probability = ' + str(prob))
            #     f.write('\n')
            #     f.write('elapsed time: ' + str(elapsed_time-start_time))
            #     f.write('\n')
            #     f.write('='*40)
            #     f.write('\n')


        bound = bound + 1
        prob = graph.model_check(model, model_name, prism, csl_prop)
        
    f.close()

f.close()
