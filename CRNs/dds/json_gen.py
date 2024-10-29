#!/usr/bin/python3

import json, os

n = 6

folder_path = "./n" + str(n)
file_path = folder_path + "/dds_" + str(n) + ".json" 
csl_path = folder_path + "/dds_" + str(n) + ".csl" 
with open(csl_path, "r") as f:
    csl_property = f.readline()

if not os.path.exists(folder_path): 
    os.makedirs(folder_path)

json_file = {}
json_file["model_path"] = "./CRNs/dds/n" + str(n) + "/dds_" + str(n) + ".sm"
json_file["jani_path"] = "./CRNs/dds/n" + str(n) + "/dds_" + str(n) + ".jani"
json_file["model_name"] = "dds_" + str(n)
json_file["csl_property"] = csl_property 
json_file["target_variables"] = ["disc_cluster0","disc_cluster1","disc_cluster2","disc_cluster3","disc_cluster4","disc_cluster5",] 
json_file["target_variables"].extend(["disc_controller0", "disc_controller1"])
json_file["target_variables"].extend(["processor"])
json_file["target_values"] = [str(n) for i in range(9)] 

try:
    with open(file_path, "w") as jf:
        json.dump(json_file, jf, indent=4, ensure_ascii=False)
    print(f"JSON data saved to '{file_path}' successfully.")
except IOError:
    print(f"Unable to save JSON data to '{file_path}'.")

