#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
from Graph import Node, Graph
from Utils import get_total_outgoing_rate, get_reaction_rate
import sys
import json
from Parser import Parser
import random
from math import log


def SSA_SS(model, N, max_time):
    ss = [0] * len(model.get_initial_state())
    for i in range(N):
        if i%100 == 0:
            print(i)
        curr_state = model.get_initial_state()

        time = 0
        while(True):
            out_rate = get_total_outgoing_rate(curr_state, model)
            rand1 = random.uniform(0, 1)
            tau = float(1.0 / float(out_rate)) * log(float(1.0 / rand1))
            time = time + tau
            if time > max_time:
                ss = [ss[i] + curr_state[i] for i in range(len(ss))]
                break
    
            rand2 = random.uniform(sys.float_info.min , out_rate)
            j = 0
            temp_sum = 0
            reaction_rates = [get_reaction_rate(curr_state, model, r_index) for r_index, _ in enumerate(model.get_reactions_vector())]
            while temp_sum <= rand2:
                temp_sum = temp_sum + reaction_rates[j]
                j = (j + 1)  
            
            j = j - 1
            next_state = [None] * len(curr_state)
            for k, s in enumerate(curr_state):
                next_state[k] = s + model.get_reactions_vector()[j][k]
            curr_state = tuple(next_state)
            
    return ([ss[i]/N for i in range(len(ss))])

def SSA(model, N, max_time, target_index, target_value):
    p_count = 0
    returned_freq_dict = {}
    returned_pop_dict = {}
    for i,r in enumerate(model.get_reactions_vector()):
            returned_freq_dict[i] = 0
    for i,s in enumerate(model.get_species_tuple()):
        returned_pop_dict[s] = 0

    for i in range(N):
        if i%100 == 0:
            print(i)
        curr_state = model.get_initial_state()
        r_count_dict = {}
        pop_dict = {}
        for i, r in enumerate(model.get_reactions_vector()):
            r_count_dict[i] = 0
        for i, s in enumerate(model.get_species_tuple()):
            pop_dict[s] = curr_state[i]
        
        time = 0
        state_count = 1
        while(True):
            out_rate = get_total_outgoing_rate(curr_state, model)
            rand1 = random.uniform(0, 1)
            tau = float(1.0 / float(out_rate)) * log(float(1.0 / rand1))
            time = time + tau
            if time > max_time:
                break
    
            rand2 = random.uniform(sys.float_info.min , out_rate)
            j = 0
            temp_sum = 0
            reaction_rates = [get_reaction_rate(curr_state, model, r_index) for r_index, _ in enumerate(model.get_reactions_vector())]
            while temp_sum <= rand2:
                temp_sum = temp_sum + reaction_rates[j]
                j = (j + 1)  
            
            j = j - 1
            r_count_dict[j] = r_count_dict[j] + 1
            next_state = [None] * len(curr_state)
            for k, s in enumerate(curr_state):
                next_state[k] = s + model.get_reactions_vector()[j][k]
            curr_state = tuple(next_state)
            state_count = state_count + 1
            for k, s in enumerate(pop_dict):
                pop_dict[s] = pop_dict[s] + curr_state[k]

            # if curr_state[target_index] >= target_value:
            #     p_count = p_count + 1
            #     break
    
        sum = 0
        for j, e in enumerate(r_count_dict):
            sum = sum + r_count_dict[j]
        if sum > 0:
            for j, e in enumerate(r_count_dict):
                r_count_dict[j] = r_count_dict[j] / sum
        else:
            print('SUM = 0 !!!')
        
        for j, e in enumerate(pop_dict):
            pop_dict[e] = pop_dict[e] / state_count

        for j, e in enumerate(pop_dict):
            returned_pop_dict[e] = returned_pop_dict[e] + pop_dict[e]/N
        
        for j, e in enumerate(r_count_dict):
            returned_freq_dict[j] = returned_freq_dict[j] + r_count_dict[j]/N

    return (returned_pop_dict, returned_freq_dict, float(p_count/N))


f = open(sys.argv[1])
json_data = json.load(f)

model_path = json_data['model_path']
target_var = json_data['target_variable']
target_value = int(json_data['target_value'])
model = Parser(model_path)
target_index = model.species_to_index_dict[target_var]
print(model.get_species_tuple())
print(SSA_SS(model, 1000, 1000))
# a, b, p = SSA(model, 2000, 50000, target_index, target_value)
# print(a)
# print(b)
# print(p)