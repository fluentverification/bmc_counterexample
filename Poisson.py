#! /home/ubu/projects/probmc/storm/pycarl/env/bin/python3
#
from scipy.stats import poisson
import math
from collections import OrderedDict
from Probability import poisson_at_least_k

print(poisson_at_least_k(model_name="enzym", poisson_step=1, division_factor=4))
# for i in range(0, 29, 3):
#     # print(poisson.cdf(k = i, mu = 300))
#     print(i)
quit()
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
                    return None
                    # raise Exception("Excpetion - poisson cdf cannot be called with k = " + str(max_k) + " and mu = " + str(r))
                    # prob = cut_off
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

def avg_discrete_probability(model_name):
    if "enzym" in model_name:
         avg_prob = [0.25042114778239605, 0.22675486210780674, 0.022675486210879, 0.2507018418572009, 0.2267696927642055, 0.022676969276518397]
    elif "motil" in model_name:
        avg_prob = [0.05235534034146383, 0.0010514715116720988, 0.3796802335433258, 0.0015577805187723377, 0.1396828812847144, 0.0012242941019240574, 0.11826681765312089, 0.013968288128471573, 0.07991553479830538, 0.06674265732859438, 0.08203439279345018, 0.06352030799621344]
    elif "yeast" in model_name:
        avg_prob = [0.0002670219293031733, 0.0007073052188598062, 0.14853409596056194, 0.01128223474863973, 0.28526512446924307, 0.16481332819257188, 0.16356762814840484, 0.22556326133242866]
    elif "circuit" in model_name:
        avg_prob = [0.1693465200442219, 0.18668013533779118, 0.00025490032325493006, 0.05070953841475635, 0.00025490032325493006, 0.04998500666371727, 0.20508683770567976, 0.18668013533779118, 0.04998500666371727, 0.05020050235260115, 0.04998500666371727, 1.1204974353254456E-05, 0.0003855232895298863, 0.00022556147754557233, 0.00020922042808897793]
    
    return avg_prob