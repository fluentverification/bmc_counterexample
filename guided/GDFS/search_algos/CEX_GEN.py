import time
from math import log, pow

def CEX_GEN(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, init_bound, init_threshold):
     start_time = time.time()
     threshold_log_step = -1
     flag = False
     bound = current_max_bound = init_bound
     flag = False
     threshold = log(init_threshold, 10)
     diag = Graph()
     while (True):
        print('-' * 50)
        print('bound = ' + str(bound))
        print('probability threshold = ' + str(pow(10, threshold)) + '\n')
        diag = BMC_DFS(model, target_var, target_index, target_value, mc_step, model_name, prism_bin, csl_prop, bound, threshold, diag)
        if signal:
            return
        if flag:
            bound = bound + 1
            flag = False
        else:
            threshold = threshold + threshold_log_step
            bound = init_bound
            flag = True