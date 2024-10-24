#!/usr/bin/python3
import sys, os

n = 2
lambda_ = 1.0 
mu = 2.0

folder_path = "./n" + str(n)
file_path = folder_path + "/dds_" + str(n) + ".sm" 

if not os.path.exists(folder_path): 
    os.makedirs(folder_path)

f = open(file_path, "w")

f.write("ctmc")
f.write("\n")
f.write("\n")

f.write("const double lambda = " + str(lambda_) + " ;")
f.write("\n")
f.write("const double mu = " + str(mu) + " ;")
f.write("\n")
f.write("\n")

f.write("module dds_" + str(n))
f.write("\n")
f.write("\n")

for i in range(6):
    f.write("disc_cluster" + str(i) + " : [0.." + str(2*n) + "] init 0;")
    f.write("\n")
f.write("\n")

f.write("disc_controller0 : [0.." + str(n) + "] init 0;")
f.write("\n")
f.write("disc_controller1 : [0.." + str(n) + "] init 0;")
f.write("\n")
f.write("\n")

f.write("processor : [0.." + str(n) + "] init 0;")
f.write("\n")
f.write("\n")

for i in range(6):
    f.write("[cluster" + str(i) + "_failing] (disc_cluster" + str(i) + " < " + str(2*n) + ") : ((" + str(2*n) + " - disc_cluster" + str(i) + ") * labmda) -> (disc_cluster" + str(i) + " = disc_cluster" + str(i) + " + 1);") 
    f.write("\n")
f.write("\n")

f.write("[controller0_failing] (disc_controller0 < " + str(n) + ")  : ((" + str(n) + " - disc_controller0) * 3 * labmda) -> (disc_controller0 = disc_controller0 + 1);") 
f.write("\n")
f.write("[controller1_failing] (disc_controller1 < " + str(n) + ")  : ((" + str(n) + " - disc_controller1) * 3 * labmda) -> (disc_controller1 = disc_controller1 + 1);") 
f.write("\n")
f.write("\n")

f.write("[processor_failing] (processor < " + str(n) + ")  : ((" + str(n) + " - processor) * 3 * labmda) -> (processor = processor + 1);") 
f.write("\n")
f.write("\n")

for i in range(6):
    f.write("[cluster" + str(i) + "_repair] (disc_cluster" + str(i) + " > 0 ) : (mu) -> (disc_cluster" + str(i) + " = disc_cluster" + str(i) + " - 1);") 
    f.write("\n")
f.write("\n")

f.write("[controller0_repair] (disc_controller0 > 0)  : (mu) -> (disc_controller0 = disc_controller0 - 1);") 
f.write("\n")
f.write("[controller1_repair] (disc_controller1 > 0)  : (mu) -> (disc_controller1 = disc_controller1 - 1);") 
f.write("\n")
f.write("\n")

f.write("[processor_repair] (processor > 0)  : (mu) -> (processor = processor - 1);") 
f.write("\n")
f.write("\n")

f.write("endmodule")

f.close()


