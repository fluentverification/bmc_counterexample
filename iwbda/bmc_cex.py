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
def bmc_dq(graph, level, step, base, property_var, max_search_bound, path_dict):
	if level == base + step: 
		bound = 1
		flag = True
		while (flag):
			solver = Solver()
			solver.add(get_initial_state(model))

			for j in range(1, bound+1):
				solver.add(get_encoding(model, j))

			x = Int(property_var + '.' + str(bound))
			property_constraint = (x==(base+step))

			curr_path_list = path_dict[(level-step, bound)]
			for path in curr_path_list:
				solver.add(exclude_path(path))
			
			if (solver.check(property_constraint) == sat): 
				flag = False
				path = solver.model()
				graph.add_path(path)
				curr_path_list = path_dict[(level-step, bound)]
				curr_path_list.append(path)
				path_dict[(level-step, bound)] = curr_path_list

			bound = bound + 1
			if bound>=max_search_bound: 
				break

	else: 
		bmc_dq(graph, level-step, step, base, property_var, max_search_bound, path_dict)

		curr_init_list = []
		for n in graph.node_list: 
			for s in n.var_dict:
				if s == property_var:
					if n.var_dict[s] == (level-step):
						curr_init_list.append(n)

		

		init_const = []
		for k, n in enumerate(curr_init_list):
			var_values = []
			for s in n.var_dict: 
				x = Int (s + '.0')
				temp_const = (x == n.var_dict[s])
				var_values.append(temp_const)
			init_const.append(And(var_values))
		init_const = Or(init_const)
			
		bound = 1
		flag = True
		while (flag):
			solver = Solver()
			solver.add(init_const)

			for j in range(1, bound+1):
				solver.add(get_encoding(model, j))

			x = Int(property_var + '.' + str(bound))
			property_constraint = (x==(level))

			curr_path_list = path_dict[(level-step, bound)]
			for path in curr_path_list:
				solver.add(exclude_path(path))

			if (solver.check(property_constraint) == sat): 
				#print('here')	
				flag = False
				path = solver.model()
				graph.add_path(path)
				curr_path_list = path_dict[(level-step, bound)]
				curr_path_list.append(path)
				path_dict[(level-step, bound)] = curr_path_list

			bound = bound+1
			if bound>max_search_bound: 
				break


##################################################################
##################################################################
##################################################################

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


prob = 0
max_search_bound = 50

path_dict = {}
for i in range(base, property_val, step):
	for j in range(1, max_search_bound):
		path_dict[(i, j)] = []

count = 0
start_time = time.time()
seed_path_total_time = 0
scaffold_total_time = 0

with open('./results/' + model_name + '/' + model_name + '.results', mode = 'w', encoding= 'ascii') as f:
	f.truncate()
	while prob <= prob_bound:
		
		seed_path_start_time = time.time()	
		bmc_dq(graph, property_val, step, base, property_var, max_search_bound, path_dict)
		seed_path_end_time = time.time()
		seed_path_total_time = seed_path_total_time + (seed_path_end_time-seed_path_start_time)
		count = count + 1
		prob = graph.model_check(model, model_name, prism, csl_prop)
		if count>2:
			count = 0
			scaffold_start_time = time.time()
			count, prob1, terminate = construct_path_3(graph, model, model_name, prism, csl_prop, cp_bound, count, prob, prob_bound, property_var, property_val, mc_step)
			scaffold_end_time = time.time()
			scaffold_total_time = scaffold_total_time + (scaffold_end_time-scaffold_start_time)
			prob = graph.model_check(model, model_name, prism, csl_prop)
		elapsed_time = time.time()
		print('# of nodes: ' + str(len(graph.node_list)))
		print('# of edges: ' + str(len(graph.edge_list)))
		print('probability = ' + str(prob))
		print('elapsed time: ' + str(elapsed_time-start_time))
		print('='*40)

		
		f.write('# of nodes: ' + str(len(graph.node_list)))
		f.write('\n')
		f.write('# of edges: ' + str(len(graph.edge_list)))
		f.write('\n')
		f.write('probability = ' + str(prob))
		f.write('\n')
		f.write('elapsed time: ' + str(elapsed_time-start_time))
		f.write('\n')
		f.write('='*40)
		f.write('\n')
		if (elapsed_time-start_time) > 1800:
			break
	f.close()
	print("seed path time: " + str(seed_path_total_time))
	print("scaffold time: " + str(scaffold_total_time)) 

f.close()


