#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#

import json
import sys, os
import Model_Bounded_prob

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


Model_Bounded_prob.CEX_GEN(json_data)

