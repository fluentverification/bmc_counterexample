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
prism = json_data['prism_binary']
prob_bound = float(json_data['probability_bound'])
property_var = json_data['property_variable']
property_val = int(json_data['property_value'])
mc_step = int(json_data['model_check_step'])
scaffold_bound_limit = int(json_data['scaffold_bound_limit'])
scaffold_count_limit = int(json_data['scaffold_count_limit'])
#base = int(json_data['base'])
#step = int(json_data['step'])
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
start_time = time.time()

with open('./results/' + model_name + '/' + model_name + '.results', mode = 'w', encoding= 'ascii') as f:
    f.truncate()
	
    bound = starting_bound
    while prob <= prob_bound:
		 
        solver = Solver()
        solver.add(get_initial_state(model))

        for j in range(1, bound+1):
            solver.add(get_encoding(model, j))
        
        #for j in range(1, bound):
        #    solver.add(loop_constraint(model, j))

        x = Int(property_var + '.' + str(bound))
        property_constraint = (x==property_val)
        #print("here")	
        while (solver.check(property_constraint) == sat):
            count = count + 1
            path = solver.model()
            solver.add(exclude_path(path))
            graph.add_path(path)
            if (count % mc_step == 0):
                prob = graph.model_check(model, model_name, prism, csl_prop)
            if (count % 3 == 0):
                scaffold(graph, model, scaffold_bound_limit, scaffold_count_limit, property_var, property_val)
            print('# of nodes: ' + str(len(graph.node_list)))
            print('# of edges: ' + str(len(graph.edge_list)))
            print('probability = ' + str(prob))
            print('='*40)
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
