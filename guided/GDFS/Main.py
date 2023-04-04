from Utils import get_min_max
from Search_Algos import ggdfs_prob, ggdfs_prob_trace_back, ggdfs_trace_back, ggdfs
from Graph_ import Node, Edge, Graph
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
init_node = Node()
init_node.initial_state = True
init_node.var_values = model.get_initial_state()
init_node.reachability_probability = 1.0
graph = Graph()
graph.add_node(init_node)
graph_size = 1
prob_thresh = 1.0
flag1 = False
flag_prob = False
bound_changed = False
min_max = get_min_max(model, bound, target_index, target_value, max_comb)
while (True):
    print('='*50)
    print("bound = " + str(bound))
    print("prob threshold = " + str(prob_thresh))
    if bound_changed:
        min_max = get_min_max(model, bound, target_index, target_value, max_comb)
    flag2 = ggdfs(graph=graph, start_node=init_node, model=model, target_index=target_index, target_value=target_value, min_max=min_max)
    #flag2 = ggdfs_trace_back(graph=graph, start_node=init_node, model=model, target_index=target_index, target_value=target_value, min_max=min_max)
    # flag2 = ggdfs_prob(graph=graph, start_node=init_node, model=model, target_index=target_index, target_value=target_value, prob_thresh=prob_thresh, min_max=min_max)
    # flag2 = ggdfs_prob_trace_back(graph=graph, start_node=init_node, model=model, target_index=target_index, target_value=target_value, prob_thresh=prob_thresh, min_max=min_max)
    #reduce the threshold until you find the first trace for the minimum bound. Afterwards
    #it alternates between reducing the threshold or increasing the bound
    # if flag1 or flag2:
    #     flag1 = True
    #     if bound_changed:
    #         prob_thresh = prob_thresh / 2.0
    #         bound_changed = False
    #     else:
    #         bound = bound + 1
    #         bound_changed = True
    # else:
    #     prob_thresh = prob_thresh / 2.0
    bound_changed = True
    bound = bound + 1
    print("state-space size = " + str(len(graph.nodes) + len(graph.edges)))
    if flag_prob or flag2:
        flag_prob = True
        print("probability = " + str(graph.check_probability(model = model, file_name_prefix= model_name, prism_bin=prism_bin, csl_prop=csl_prop)))
    else:
        print("probability = 0")
    current_time = time.time()
    print("time = " + str(current_time - start_time))
    print('='*50)
    graph = Graph()
    graph.add_node(init_node)
    if (current_time - start_time > 3600):
        break