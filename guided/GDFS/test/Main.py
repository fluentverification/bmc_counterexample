# from Utils import get_constraints_range
# from Search_Algos import ggdfs_prob, ggdfs, post_process, XBF, Guided_Greedy_DFS, Guided_Greedy_DFS_probability
# from search_algos.CEX_GEN import CEX_GEN
# from Graph import Node, Edge, Graph
# from Model import Model
import time
import json
import sys, os
from Graph.Graph import Graph
from Graph.utils import check_probability
from Graph.Node import Node
from Graph.Edge import Edge
from Search_Algos.BMC_GDFS.CEX_GEN import CEX_GEN
from Parser import Parser




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

model = Parser(model_path)
# diag = Graph()
# node = Node()
# node.var_values = model.get_initial_state()
# node2 = Node()
# node2.var_values = (0, 49, 1, 1, 50, 0)
# edge = Edge()
# edge.src = node
# edge.dst = node2
# edge.reaction = 0
# diag.add_edge(edge)

# print (enabled(model, diag, node, 0, -1000))
# print('=' * 50)
# print(model.get_reactions_vector())
target_index = model.species_to_index_dict[target_var]

if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)

# start_time = time.time()

# print("="*50 + "\n")
# print("Starting... \n")


CEX_GEN(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, 1E-20)