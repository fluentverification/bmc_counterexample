#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3


import itertools
from z3 import Solver, Int, And, Or, sat, simplify, ToInt
from Parser import Parser
import time
import subprocess
from Min_Max import get_min_max_avg_prob
import math
from JANI_parser import JSON_Parser


#
def CEX_GEN(json_data):
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
    
    #generate subsets (keys for the min_max dictionary)
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    # max_comb = 1
    subsets = get_subsets(index_tuple, 1, max_comb)
    print("\nStarting...")
    
    #
    min_max_dict = {}
    for s in subsets:
        temp = 0
        for e in s:
            temp = temp + model.get_initial_state()[e]
        min_max_dict[s] = [temp, temp]
    
    ########### parameters ###########
    #
    thresh = -1
    division_factor = 100
    engine = "automatic"
    #
    N = math.floor(math.log(K))

    while (True):
        print("thresh = " + str(thresh))
        before = time.time()
        
        flag, length, min_max_dict = get_min_max_avg_prob(model, thresh, target_index, target_value, subsets, min_max_dict, model_name, division_factor, N)
        print("Generating min_max dictonary took " + str(time.time() - before) + "seconds.")
        #
        if flag:
            N = math.floor(math.log(length))
            JSON_Parser(model, model_name, K, jani_path, min_max_dict)
            print("Calling Storm to calculate the probability... \n\n")
            
            #running storm on the produced output
            stdout_result = subprocess.run([storm_bin, "--jani", "./results/" + model_name + "/bounds/" + model_name + "_" + str(K) + ".jani",'--engine', engine, '--prop', csl_prop_lb], stdout=subprocess.PIPE)
            stdout_result = stdout_result.stdout.decode('utf-8')
            print(stdout_result)
        else:
            print("No additional witnesses found for threshold " + str(thresh))
        #
        
        thresh = thresh - 1
        K = K + 1
        print("\n \nRunning time: " + str(time.time() - start_time) + " seconds")
        print("=" * 50)
#


#returns all the subsets of a tuple with cardinality in range [L, U]
#returns a list of tuples where each tuple is a subset
def get_subsets (tuple, L, U):
    if U>len(tuple):
        U = len(tuple)
    return list(itertools.chain.from_iterable(itertools.combinations(tuple, r) for r in range(L, U+1)))
#


