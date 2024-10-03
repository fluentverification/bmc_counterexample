#! /home/ubu/projects/probmc/storm-1.9/env/bin/python3
#

import json
import sys, os
import Model_Bounded_prob

#read json elements
f = open(sys.argv[1])
json_data = json.load(f)
model_name = json_data['model_name']

print("hello")
quit()

if not os.path.exists('./results'): 
    os.makedirs('./results')
if not os.path.exists('./results/' + model_name):
    os.makedirs('./results/' + model_name)
if not os.path.exists('./results/' + model_name + '/bounds'):
    os.makedirs('./results/' + model_name + '/bounds')


Model_Bounded_prob.CEX_GEN(json_data)



#collecting all the variables in the program{
variables = []
for v in prism_program.global_integer_variables:
    variables.append(v)
for v in prism_program.global_boolean_variables:
    variables.append(v)
for m in prism_program.modules:
    for v in m.integer_variables:
        variables.append(v)
    for v in m.boolean_variables:
        variables.append(v)
for v in variables:
    print(v.name)
    print(str(type(v)) + " : "  + str(v.initial_value_expression))
#}
print("then:")
#getting the initial state of the program as a vector of bools and ints{
initial_state_stormpy = [None] * len(variables)
for i, v in enumerate(variables):
    initial_state_stormpy[i] = v.initial_value_expression
if None in initial_state_stormpy:
    raise Exception("initial state was not retrieved")
initial_state = [None] * len(variables)
for i, v in enumerate(initial_state_stormpy):
    if v.has_integer_type():
        initial_state[i] = v.evaluate_as_int()
    elif v.has_boolean_type():
        initial_state[i] = v.evaluate_as_bool()
    else:
        raise Exception("one of the variables in the model has a type other than bool/int")
print(initial_state)
#}
