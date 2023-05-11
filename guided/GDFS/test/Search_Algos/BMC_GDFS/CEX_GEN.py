import time
from math import log, pow
from Graph.Graph import Graph
from Graph.utils import check_probability
from Search_Algos.BMC_GDFS.BMC_GDFS import BMC_GDFS

def CEX_GEN(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound, init_threshold):
    start_time = time.time()
    step = -1
    K = K_max = init_bound
    flag = False
    P = log(init_threshold, 10)
    diag = Graph()
    diag_size = len(diag.nodes) + len(diag.edges)
    
    while (True):
        diag = BMC_GDFS(model, K, P, diag, target_var, target_index, target_value)
        # print(P)
        # print(K)
        # print('-'*50)
        if (len(diag.nodes) + len(diag.edges) - diag_size >= mc_step):
            result = check_probability(diag, model, model_name, prism_bin, csl_prop)
            print('probability= ' + str(result))
            diag_size = len(diag.edges) + len(diag.nodes)
            print('size= ' + str(diag_size))
            print('time= ' + str(time.time()-start_time))
            print('='*30)
            if time.time()-start_time > 3600:
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