from Utils import get_constraints_range
from Search_Algos import ggdfs_prob, ggdfs, post_process, XBF, Guided_Greedy_DFS, Guided_Greedy_DFS_probability
from search_algos.BMC_DFS import BMC_DFS, trace_as_list
from Graph import Node, Edge, Graph
from Model import Model
import time
import json
import sys, os


#read json elements
f = open(sys.argv[1])
json_data = json.load(f)
model_path = json_data['model_path']
model_name = json_data['model_name']
bound = int(json_data['starting_bound'])
csl_prop = json_data['csl_property']
prism_bin = json_data['prism_binary']
target_var = json_data['target_variable']
target_value = int(json_data['target_value'])
mc_step = int(json_data['model_check_step'])
max_comb = int(json_data['max_combination'])

model = Model(model_path)
target_index = model.species_to_index_dict[target_var]

if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)

start_time = time.time()

print("="*50 + "\n")
print("Starting... \n")

###############################################################################################################
# circuit_range_list = [(3,40)]
# motil_range_list = [(10, 15), (15, 20), (20, 25), (25,30), (30,35)]
# yeast_range_list = [(100,102)]
# bnimply_range_list = [(1,5), (5,10), (10,15), (15,20), (20,25)]
# range_list = bnimply_range_list
# min_max_list = get_constraints_range(model, range_list, target_index, target_value, max_comb)
# print("\nmin_max_list generated in " + str(time.time()-start_time) + " seconds.")
###############################################################################################################
min_max_list = [[0],[0]]
print("\n--------------------")

enzym_prob_thresh = 1E-25
yeast_prob_thresh = 1E-150
motil_prob_thresh = 1E-20
circuit_rpob_thresh = 1E-20
################################################################################################
# XBF(model, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, min_max_list[0])
# Guided_Greedy_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound)
# Guided_Greedy_DFS_probability(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, 1E-20)
BMC_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, 1E-20)
################################################################################################
