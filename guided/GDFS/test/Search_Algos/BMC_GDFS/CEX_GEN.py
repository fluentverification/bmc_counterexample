import time
from math import log, pow
from Graph.Graph import Graph
from Graph.Edge import Edge
from Graph.utils import check_probability
from Search_Algos.BMC_GDFS.BMC_GDFS import BMC_GDFS
from Utils.misc import get_reaction_rate

def CEX_GEN(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound, init_threshold):
    start_time = time.time()
    step = -1
    K = K_max = init_bound
    flag = False
    P = log(init_threshold, 10)
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)
    E_count = sat_count = unsat_count = 0
    sat_time = unsat_time = 0
    nodes = {}
    func_count = 0
    total_diag_size = 0
    
    while (True):
        total_diag_size = total_diag_size + len(diag.nodes) + len(diag.edges)
        func_count = func_count + 1

        diag, E_count, sat_count, unsat_count, sat_time_temp, unsat_time_temp = BMC_GDFS(model, K, P, diag, target_var, target_index, target_value, E_count, sat_count, unsat_count, nodes)
        sat_time = sat_time + sat_time_temp
        unsat_time = unsat_time + unsat_time_temp
        
        for n in diag.nodes:
            for n_prime in diag.nodes:
                for i, r in enumerate(model.get_reactions_vector()):
                    dst_var_values = [None] * len(model.get_species_tuple())
                    for j, coefficient in enumerate(r):
                        dst_var_values[j] = n[j] + coefficient
                    dst_var_values = tuple(dst_var_values)
                    if dst_var_values == n_prime:
                        edge = Edge()
                        edge.src = diag.nodes[n]
                        edge.dst = diag.nodes[n_prime]
                        edge.reaction = i
                        edge.rate = get_reaction_rate(n, model, i)
                        diag.add_edge(edge)
        # print(P)
        # print(K)
        # print('-'*50)
        if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
            result = check_probability(diag, model, model_name, prism_bin, csl_prop)
            print('probability= ' + str(result))
            diag_size = len(diag.edges) + len(diag.nodes)
            print('size= ' + str(diag_size))
            print('time= ' + str(time.time()-start_time))
            print('-' * 10)
            print('frac_sat_time= ' + str(float(sat_time/ (time.time()-start_time))))
            print('frac_unsat_time= ' + str(float(unsat_time/ (time.time()-start_time))))
            print('E_count = ' + str(E_count) + ' | sat_count = ' + str(sat_count) + ' | unsat_count = ' + str(unsat_count))
            print('#BMC_GDFS calls = ' + str(func_count))
            print('total diag size passed = ' + str(total_diag_size))
            print('='*50)
            if time.time()-start_time > 1800:
                quit()
        
        if flag:
            K = K + 1
            if K > K_max:
                K_max = K
                flag = False
        else:
            P = P + step
            K = init_bound
            flag = True