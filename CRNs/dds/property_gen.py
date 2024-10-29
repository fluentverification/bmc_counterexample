#!/usr/bin/python3
import sys, os
import math

n = 6

folder_path = "./n" + str(n)
file_path = folder_path + "/dds_" + str(n) + ".csl" 

if not os.path.exists(folder_path): 
    os.makedirs(folder_path)

f = open(file_path, "w")

failure = ""
count = 0 
for i in range(6):
    if count != 0:
        failure = failure + " | " 
    failure = failure + "(disc_cluster" + str(i) + " >= " + str(n) + ")"
    count = count + 1

failure = failure + " | (disc_controller0 = " + str(n) + ")"
failure = failure + " | (disc_controller1 = " + str(n) + ")"

failure = failure + " | (processor = " + str(n) + ")"


prop = "P=? [ true U<=840 (" + str(failure) + ")]"
f.write(prop)
f.close()

