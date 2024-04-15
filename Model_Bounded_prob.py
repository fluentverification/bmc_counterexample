#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3


import itertools
from z3 import Solver, Int, And, Or, sat, simplify, ToInt
from Parser import Parser
import time
import subprocess
from Min_Max_ import get_min_max_species
import math
from JSON_Parser import JSON_Parser
import stormpy


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
    thresh = 0
    delta = 1
    division_factor = 1000
    engine = "automatic"
    # max_comb_species = len(model.get_species_tuple())
    max_comb_species = 1
    steps = 1
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
    min_max_dict_species = {}
    for s in subsets_species:
        temp = 0
        for e in s:
            temp = temp + model.get_initial_state()[e]
        min_max_dict_species[s] = [temp, temp]
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
                                                          min_max_prev=min_max_dict_species, 
                                                          target_index=target_index, 
                                                          target_value=target_value, 
                                                          lower_bound = lower_bound)
        print("Generating min_max dictonary took " + str(time.time() - before) + "seconds.")
        #
        
        if flag:
            # print(min_max_dict_species)
            # new_steps = math.floor(math.log(max_len))
            # if new_steps > steps:
            #     steps = new_steps
            
            dest_path = "./results/" + model_name + ".jani"
            JSON_Parser(model=model, 
                        dest_path=dest_path,  
                        jani_path=jani_path,
                        min_max_dict=min_max_dict_species)

            print("Calling Storm to calculate the probability... \n\n")
            
            #running storm on the produced output
            jani, _ = stormpy.parse_jani_model(dest_path)
            jani_property = stormpy.parse_properties_for_jani_model(csl_prop_lb, jani)
            sparse_model = stormpy.build_sparse_model(jani, jani_property)
            ss_size = int(sparse_model.nr_states) + int(sparse_model.nr_transitions)
            print("state-space size = " + str(ss_size))
            try:
                result = stormpy.check_model_sparse(sparse_model, jani_property[0], only_initial_states=True)
            except:
                print("model-checking for this iteration could not be done")
                print("=" * 50)
                thresh = thresh + delta
                continue
            filter = stormpy.create_filter_initial_states_sparse(sparse_model)
            result.filter(filter)
            print("Lower-Bound Probability = " + str(float(result.max)))
            #
        else:
            print("No additional witnesses found for threshold " + str(thresh))
        #
        
        thresh = thresh + delta

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


