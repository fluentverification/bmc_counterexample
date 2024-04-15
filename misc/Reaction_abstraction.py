#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#

import json
import sys, os
import time
import subprocess
from Parser import Parser
from JANI_parser import JSON_Parser_reaction


#read json elements
f = open(sys.argv[1])
json_data = json.load(f)
model_name = json_data['model_name']


if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)
if not os.path.exists('./results/' + model_name + '/bounds'):
    os.makedirs('./results/' + model_name + '/bounds')


start_time = time.time()
#parsing parameters from json elements
jani_path = json_data['jani_path']
model_name = json_data['model_name']
model_path = json_data['model_path']
K = int(json_data['starting_bound'])
csl_prop_lb = json_data['csl_property']
csl_prop_ub = json_data['csl_property_ub']
target_var = json_data['target_variable']
target_value = int(json_data['target_value'])
max_comb = int(json_data['max_combination'])
storm_bin = json_data["storm_binary"]
#
#parse the model into Parser() object
model = Parser(model_path)
target_index = model.species_to_index_dict[target_var]
#
print(model.get_reactions_vector())
print(model.get_species_tuple())
min_max_reactions = {0 : 40,
                     1 : 40,
                     2 : 40,
                     3 : 40,
                     4 : 40,
                     5 : 40
                     }


JSON_Parser_reaction(model, model_name, 0, jani_path, min_max_reactions, model.index_to_reaction_dict)
stdout_result = subprocess.run([storm_bin, "--jani", "./results/" + model_name + "/bounds/" + model_name + "_" + str(0) + ".jani",'--engine', 'automatic', '--prop', csl_prop_lb, '--exportbuild', "./test.dot"], stdout=subprocess.PIPE)
stdout_result = stdout_result.stdout.decode('utf-8')
print(stdout_result)