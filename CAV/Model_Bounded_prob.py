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
    thresh = 0
    delta = 2
    division_factor = 100
    engine = "sparse"
    max_comb_species = len(model.get_species_tuple())
    # max_comb_species = 2
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
    # for s in subsets_species:
    #     temp = 0
    #     for e in s:
    #         temp = temp + model.get_initial_state()[e]
    #     min_max_dict_species[s] = [temp, temp]
    #
    for i in range(steps):
        min_max_bound = {}
        for s in subsets_species:
            min_max_bound[s] = [None, None]
        min_max_dict_species[i] = min_max_bound

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
        print("=" * 50)

        if flag:
            # print(min_max_dict_species)
            
            JANI_Parser(model=model, 
                        model_name=model_name, 
                        jani_model=jani_model, 
                        min_max_dict=min_max_dict_species)

            # quit()
            print("Calling Storm to calculate the probability... \n\n")
            #Calculating P(A|B)
            jani, _ = stormpy.parse_jani_model("./tmp/" + model_name +  "_sink" + ".jani")
            jani_property = stormpy.parse_properties_for_jani_model(prop_lb, jani)
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
            A_B = float(result.min)
            #
            print("lower-bound = " + str(A_B))
            #Calculating P(B)
            # jani, _ = stormpy.parse_jani_model("./tmp/" + model_name +  "_sink" + ".jani")
            # jani_property = stormpy.parse_properties_for_jani_model(prop_ub, jani)
            # sparse_model = stormpy.build_sparse_model(jani, jani_property)
            # ss_size = int(sparse_model.nr_states) + int(sparse_model.nr_transitions)
            # print("state-space size = " + str(ss_size))
            # try:
            #     result = stormpy.check_model_sparse(sparse_model, jani_property[0], only_initial_states=True)
            # except:
            #     print("model-checking for this iteration could not be done")
            #     print("=" * 50)
            #     thresh = thresh + 1
            #     continue
            # filter = stormpy.create_filter_initial_states_sparse(sparse_model)
            # result.filter(filter)
            # B = 1.0 - float(result.min)
            # #
            
            # print("B = " + str(B) )
            # prob_lb = A_B * B

            # print("lower-bound probability = " + str(prob_lb))
            # prob_ub = (A_B * B) / ( (A_B * B) / ( (A_B * B) + (1 - B) ))
            # print("upper-bound probability = " + str(prob_ub))

            new_steps = math.floor(math.log(max_len))
            if new_steps > steps:
                for i in range(steps, new_steps):
                    min_max_bound = {}
                    for s in subsets_species:
                        min_max_bound[s] = [None, None]
                    min_max_dict_species[i] = min_max_bound
                steps = new_steps
        else:
            print("No additional witnesses found for threshold " + str(thresh))
        
       
        
        thresh = thresh + delta
        print("\n \nRunning time: " + str(time.time() - start_time) + " seconds")
        print("=" * 50)
#



