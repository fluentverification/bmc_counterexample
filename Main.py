#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3

import json
import sys, os
from Search_Algorithms import Bounded_DFS_Prob, Bounded_DFS, XBF, Bounded_DFS_r2
import Model_Bounded
import Model_Bounded_prob
from Parser import Parser

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


# Bounded_DFS.CEX_GEN(json_data)
# Bounded_DFS_Prob.CEX_GEN(json_data)
# XBF.CEX_GEN(json_data)
Model_Bounded_prob.CEX_GEN(json_data)
# Bounded_DFS_r2.CEX_GEN(json_data)

