#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3

import time
from Min_Max_ import get_min_max_species
import math
from JANI_Parser import JANI_Parser
import stormpy
from Utils import get_subsets


#
def CEX_GEN(model, model_name, prop_lb, prop_ub, target_variable, target_value, jani_model):
    start_time = time.time()

    ##############################
    ########## parameters ########
    #
    thresh = -1
    division_factor = 1000
    engine = "sparse"
    # max_comb_species = len(model.get_species_tuple())
    max_comb_species = 1
    steps = 2
    lower_bound = True
    poisson_step = 10
    #
    ##############################
    ##############################
    
    target_index = model.species_to_index_dict[target_variable]

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
            print(min_max_dict_species)
            new_steps = math.floor(math.log(max_len))
            if new_steps > steps:
                steps = new_steps
            
            JANI_Parser(model=model, 
                        model_name=model_name, 
                        file_suffix=0-thresh, 
                        jani_model=jani_model, 
                        min_max_dict=min_max_dict_species)

            # quit()
            print("Calling Storm to calculate the probability... \n\n")
            jani, _ = stormpy.parse_jani_model("./tmp/" + model_name +  str(0-thresh) + ".jani")
            jani_property = stormpy.parse_properties_for_jani_model(prop_lb, jani)
            sparse_model = stormpy.build_sparse_model(jani, jani_property)
            result = stormpy.check_model_sparse(sparse_model, jani_property[0], only_initial_states=True)
            filter = stormpy.create_filter_initial_states_sparse(sparse_model)
            result.filter(filter)
            print("probability = " + str(float(result.min)))
        else:
            print("No additional witnesses found for threshold " + str(thresh))
        
        thresh = thresh - 1
        print("\n \nRunning time: " + str(time.time() - start_time) + " seconds")
        print("=" * 50)
#



