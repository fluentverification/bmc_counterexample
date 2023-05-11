import time
from math import log, pow
from search_algos.BMC_DFS import BMC_DFS
from search_algos.enabled import enabled
from Graph import Graph

def CEX_GEN(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound, init_threshold):
    start_time = time.time()
    threshold_log_step = -1
    K = K_max = init_bound
    flag = False
    P = log(init_threshold, 10)
    diag = Graph()
    diag_size = 0
    while (True):
        # print('-' * 50)
        # print('bound = ' + str(bound))
        # print('probability threshold = ' + str(pow(10, threshold)) + '\n')
        diag = BMC_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, threshold, start_time, diag, diag_size)
        if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
            result = diag.check_probability(model, model_name, prism_bin, csl_prop)
            print('probability= ' + str(result))
            diag_size = len(diag.edges) + len(diag.nodes)
            print('size= ' + str(diag_size))
            print('time= ' + str(time.time()-start_time))
            print('='*30)
            if time.time()-start_time > 3600:
                quit()
        # if signal:
        #     return
        if flag:
            bound = bound + 1
            if bound > current_max_bound:
                current_max_bound = bound
                flag = False
        else:
            threshold = threshold + threshold_log_step
            bound = init_bound
            flag = True