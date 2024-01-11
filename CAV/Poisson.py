#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#
from scipy.stats import poisson
import math
from collections import OrderedDict

def poisson_cdf(model_name, cut_off, division_factor):
    
    if "enzym" in model_name:
        avg_rate = [107.199, 96.646, 9.569, 107.373, 96.796, 9.588]
    elif "motil" in model_name:
        avg_rate = [1.029, 0.022, 7.514, 0.019, 2.76, 0.015, 2.263, 0.27, 1.574, 1.23, 1.533, 1.203]
    elif "yeast" in model_name:
        avg_rate = [0.069, 0.194, 42.009, 3.259, 81.074, 46.68, 46.68, 64.428]
    elif "circuit" in model_name:
        avg_rate = [169.015, 186.398, 0.263, 50.382, 0.261, 49.758, 203.845, 185.522, 49.74, 50.089, 49.564, 0.013, 0.378, 0.236, 0.215]
    elif "single" in model_name:
        avg_rate = [1.0, 0.025]
        
    max_k_vector = [None] * len(avg_rate)
    for i, r in enumerate(avg_rate):
        max_k = 0
        prob = 0
        if r == 0:
            max_k_vector[i] = -1
        else:
            while (prob > cut_off):
                try:
                    prob = math.log(1-poisson.cdf(k = max_k, mu=r), division_factor)
                except:
                    raise Exception("Excpetion - poisson cdf cannot be called with k = " + str(max_k) + " and mu = " + str(r))
                    prob = cut_off
                max_k = max_k + 1
            max_k_vector[i] = max_k - 1
    
    poisson_ = []
    for i, e in enumerate(max_k_vector):
        dict = OrderedDict()
        if e == -1:
            dict = None
        else:
            for ii in range(e+1):
                dict[ii] = poisson.cdf(k = ii, mu=avg_rate[i])     
        poisson_.append(dict)
    
    return poisson_


def poisson_cdf_step(model_name, step, cut_off, division_factor):
    if step == 1:
        raise Exception("step provide to poisson_cdf is 1. Value greater than 1 expected.")
    
    if "enzym" in model_name:
        avg_rate = [107.199, 96.646, 9.569, 107.373, 96.796, 9.588]
    elif "motil" in model_name:
        avg_rate = [1.029, 0.022, 7.514, 0.019, 2.76, 0.015, 2.263, 0.27, 1.574, 1.23, 1.533, 1.203]
    elif "yeast" in model_name:
        avg_rate = [0.069, 0.194, 42.009, 3.259, 81.074, 46.68, 46.68, 64.428]
    elif "circuit" in model_name:
        avg_rate = [169.015, 186.398, 0.263, 50.382, 0.261, 49.758, 203.845, 185.522, 49.74, 50.089, 49.564, 0.013, 0.378, 0.236, 0.215]

    max_k_vector = [None] * len(avg_rate)
    for i, r in enumerate(avg_rate):
        max_k = 0
        prob = 0
        if r == 0:
            max_k_vector[i] = -1
        else:
            while (prob > cut_off):
                prob = math.log(1-poisson.cdf(k = max_k, mu=r), division_factor)
                max_k = max_k + step
            max_k_vector[i] = max_k - step
    
    poisson_approximation = []
    for i, e in enumerate(max_k_vector):
        dict = OrderedDict()
        if e == -1:
            dict = None
        else:
            l = 0
            r = step - 1
            for ii in range (0, e, step):
                tuple_ = (l, r)
                sum = 0
                for iii in range(l, r+1):
                    sum = sum + poisson.cdf(k = iii, mu=avg_rate[i])
                dict[tuple_] = math.log((sum / step), division_factor)
                l = l + step
                r = r + step
        poisson_approximation.append(dict)
    
    return poisson_approximation
