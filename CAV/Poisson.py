#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#
from scipy.stats import poisson
import math


def poisson_cdf(model_name, step, cut_off, division_factor):
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
        while (prob > cut_off):
            prob = math.log(1-poisson.cdf(k = max_k, mu=avg_rate[i]), division_factor)
            max_k = max_k + step
        max_k_vector[i] = max_k - step
    
    dict = {}
    l = 0
    r = step - 1
    for i, e in enumerate(max_k_vector):
        tuple_ = (l, r)
        sum = 0
        for ii in range(l, r+1):
            sum = sum + poisson.cdf(k = max_k, mu=avg_rate[i])
        dict[tuple_] = math.log((sum / step), division_factor)
        l = l + step
        r = r + step
    
    return dict


    
    return max_k_vector
    

print(poisson_cdf("enzym", 5, -2, 10))
