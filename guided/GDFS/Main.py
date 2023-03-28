from Utils import get_min_max
from Search_Algos import gdfs_prob
from Graph_ import Node, Edge, Graph
from Models import test, enzymatic_futile_cycle, motility_regulation, yeast_polarization
import time

prism_bin = "../../../prism-4.7-src/prism/bin/prism"
#sys.setrecursionlimit(1000000)
#Enzymatic Futile cycle
# csl_prop = "P=? [true U<=100 S5=40]"
# target_index = 4
# target_value = 40
# #original bound is 19
# bound = 19
# model = enzymatic_futile_cycle()
# prob_thresh = 1E-29
# prob_thresh = 0

#Yeast Polarization

# csl_prop = "P=? [true U<=20 Gbg=50]"
# target_index = 5
# target_value = 50
# #original bound is 100
# bound = 100
# model = yeast_polarization()
# prob_thresh = 1.0

#Motility Regulation
csl_prop = "P=? [true U<=10 CodY=19]"
target_index = 6
target_value = 19
#original bound is 9
bound = 9
model = motility_regulation()
prob_thresh = 1.0

start_time = time.time()
init_node = Node()
init_node.initial_state = True
init_node.var_values = model.get_initial_state()
init_node.bound = bound
init_node.reachability_probability = 1.0
graph = Graph()
graph.add_node(init_node)
#gdfs(graph=graph, start_node=init_node, model=model, prism_bin=prism_bin, csl_prop=csl_prop, target_index=target_index, target_value=target_value)

f1 = False
b_flag = False
even_odd = 0
min_max = get_min_max(model, bound, target_index, target_value, len(model.get_initial_state()))
while (True):
    even_odd = even_odd + 1
    if b_flag:
        min_max = get_min_max(model, bound, target_index, target_value, len(model.get_initial_state()))
    f2 = gdfs_prob(graph=graph, start_node=init_node, model=model, prism_bin=prism_bin, csl_prop=csl_prop, target_index=target_index, target_value=target_value, prob_thresh=prob_thresh, start_time=start_time, min_max=min_max)
    if f1 or f2:
        f1 = True
        if even_odd%2 == 0:
            prob_thresh = prob_thresh / 2.0
            b_flag = False
        else:
            bound = bound + 1
            b_flag = True
    else:
        prob_thresh = prob_thresh / 2.0
    print(prob_thresh)
    print(graph.check_probability(model, 'test', prism_bin, csl_prop))
    current_time = time.time()
    if (current_time - start_time > 440000):
        break

