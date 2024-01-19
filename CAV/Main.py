#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#

import json
import sys, os, subprocess
import Model_Bounded_prob
from PRISM_Parser import PRISM_Parser

#read json elements
f = open(sys.argv[1])
json_data = json.load(f)
#

#tmp directory is used to hold the jani model converted from prism file and the generated jani models
if not os.path.exists("./tmp"):
    os.makedirs("./tmp")
else: 
    for filename in os.listdir('./tmp'):
        if os.path.isfile(os.path.join('./tmp', filename)):
            os.remove(os.path.join('./tmp', filename))
#

#Parsing the PRISM model 
model = PRISM_Parser(json_data["model_path"])
#

#running storm-conv to get the Jani model
#./storm-conv --prism ./enzym.sm --prop "P=? [true U<=100 s5=40]" --tojani enzym.jani -pc
jani_path = "./tmp/" + str(json_data["model_name"]) + ".jani"
subprocess.run([json_data["storm-conv"], 
                "--prism", json_data["model_path"],
                "--prop", json_data["csl_property_lb"],
                "--tojani", jani_path,
                "-pc"
                ], 
                stdout=subprocess.PIPE)
#

#Running the algorithm
csl_property_lb = json_data["csl_property_lb"]
csl_property_ub = json_data["csl_property_ub"]
target_variable = json_data["target_variable"]
target_value = json_data["target_value"]
Model_Bounded_prob.CEX_GEN(model = model, 
                           model_name = json_data["model_name"],
                           prop_lb = csl_property_lb,
                           prop_ub = csl_property_ub,
                           target_variable = target_variable,
                           target_value = target_value,
                           jani_model = jani_path)
#



#deleting the tmp directory after the program finishes
for filename in os.listdir('./tmp'):
    if os.path.isfile(os.path.join('./tmp', filename)):
        os.remove(os.path.join('./tmp', filename))
os.rmdir("./tmp")