#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3


import itertools
from z3 import Solver, Int, And, Or, sat, simplify, ToInt
from Parser import Parser
import time
import subprocess
from Min_Max_ import get_min_max_species
import math
from JANI_parser import JSON_Parser_bound, JSON_Parser


#
def CEX_GEN(json_data):
    
    start_time = time.time()
    #parsing parameters from json elements
    jani_path = json_data['jani_path']
    model_name = json_data['model_name']
    model_path = json_data['model_path']
    csl_prop_lb = json_data['csl_property']
    csl_prop_ub = json_data['csl_property_ub']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    storm_bin = json_data["storm_binary"]
    #
    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    #

    ##############################
    ########## parameters ########
    #
    thresh = -1
    division_factor = 100
    engine = "sparse"
    max_comb_species = len(model.get_species_tuple())
    max_comb_species = 2
    steps = 2
    lower_bound = True
    poisson_step = 10
    #
    ##############################
    ##############################

    #generate subsets (keys for the min_max dictionary)
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    subsets_species = get_subsets(index_tuple, 1, max_comb_species)
    print("\nStarting...")
    #
     
    #initializing the min_max for each subset based on the population of 
    #species in the initial state of the model
    min_max_species = {}
    # for i in range(steps):
    #     min_max_bound = {}
    #     for s in subsets_species:
    #         temp = 0
    #         for e in s:
    #             temp = temp + model.get_initial_state()[e]
    #         min_max_bound[s] = [temp, temp]
    #     min_max_species[i] = min_max_bound
    for s in subsets_species:
        temp = 0
        for e in s:
            temp = temp + model.get_initial_state()[e]
        min_max_species[s] = [temp, temp]
    #

    while (True):
        print("thresh = " + str(thresh))
        before = time.time()
        flag, max_len, min_max_dict_species = get_min_max_species(model=model, 
                                                          model_name=model_name, 
                                                          prob_thresh=thresh, 
                                                          num_steps = steps, 
                                                          division_factor=division_factor,
                                                          poisson_step= poisson_step, 
                                                          subsets=subsets_species, 
                                                          min_max_prev=min_max_species, 
                                                          target_index=target_index, 
                                                          target_value=target_value, 
                                                          lower_bound = lower_bound)
        print("Generating min_max dictonary took " + str(time.time() - before) + "seconds.")
        #
        
        if flag:
            new_steps = math.floor(math.log(max_len))
            if new_steps > steps:
                steps = new_steps
            # if new_steps > steps:
            #     for i in range(steps, new_steps):
            #         min_max_bound = {}
            #         for s in subsets_species:
            #             temp = 0
            #             for e in s:
            #                 temp = temp + model.get_initial_state()[e]
            #             min_max_bound[s] = [temp, temp]
            #         min_max_species[i] = min_max_bound
            
            JSON_Parser(model=model, 
                        model_name=model_name, 
                        K=(0-thresh), 
                        jani_path=jani_path, 
                        min_max_dict=min_max_dict_species)
                        # max_len = max_len)
            print("Calling Storm to calculate the probability... \n\n")
            
            #running storm on the produced output
            if lower_bound:
                stdout_result = subprocess.run([storm_bin, "--jani", "./results/" + model_name + "/bounds/" + model_name + "_" + str(0 - thresh) + ".jani",'--engine', engine, '--prop', csl_prop_lb], stdout=subprocess.PIPE)
            else:
                stdout_result = subprocess.run([storm_bin, "--jani", "./results/" + model_name + "/bounds/" + model_name + "_" + str(0 - thresh) + ".jani",'--engine', engine, '--prop', csl_prop_ub], stdout=subprocess.PIPE)
            stdout_result = stdout_result.stdout.decode('utf-8')
            print(stdout_result)
            #
        else:
            print("No additional witnesses found for threshold " + str(thresh))
        #
        
        thresh = thresh - 1
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


