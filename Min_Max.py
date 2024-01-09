from z3 import Solver, Int, And, Or, sat
import math

def get_min_max_avg_prob(model, prob_thresh, target_index, target_value, subsets, min_max_prev, model_name, division_factor, N):

    min_max = {}
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    vars = []
    for i in range(N + 1):
        for ii, r in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            vars.append(x)
        
    constraints = []
    
    #first constraint (n_0_0,n_0_1,n_1_0,n_1_1,... >=0)
    for i in vars:
        constraints.append(i>=0)
    #
   
    #second constraint (probability threshold)
    prob = 0
    avg_prob_ = avg_prob(model_name)
    for i, _ in enumerate(model.get_reactions_vector()):
        for ii in range(N + 1):
            x = Int("n_" + str(ii) + "_" + str(i))
            prob = prob + (math.log(avg_prob_[i], division_factor) * x)

    constraints.append((prob>prob_thresh))
    #
    
    #third constraint (reaching the target)
    vars_ = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            for ii in range(N+1):
                vars_.append([Int("n_" + str(ii) + "_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars_:
        sum = sum + i[0]*i[1]

    constraints.append(sum==target_value)
    #

    #fourth constraint (species population>=0)
    species_pops = [None] * (N + 1)
    initial_valuation = list(model.get_initial_state())
    species_pops[0] = initial_valuation
    for i in range(1, N+1):
        valuation = [None] * len(model.get_species_tuple())
        for ii, s in enumerate(model.get_species_tuple()):
            valuation[ii] = initial_valuation[ii]
            for iii, r in enumerate(model.get_reactions_vector()):
                x = Int("n_" + str(i-1) + "_" + str(iii))
                if r[ii] > 0:
                    valuation[ii] = valuation[ii] + (x * r[ii])
        species_pops[i] = valuation
        initial_valuation = valuation
    for i, v in enumerate(species_pops):
        if i > 0:
            for ii in v:
                constraints.append(ii>=0)
    #

    #
    
 
    solver = Solver()
    solver.add(And(constraints))

    #upper-bound
    length = 0
    flag_u = False
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars__ = [None] * (N + 1)
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars__[state_index] == None:
                        vars__[state_index] = []
                        vars__[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars__[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            max_value = 0
            for ii in e:
                max_value = max_value + model.get_initial_state()[ii]
            for ii in range(N + 1):
                for iii in vars__[ii]:
                    for iv in e:
                        max_value = max_value + (model.get_reactions_vector()[iii[0]][iv] * iii[1])
                if max_value > min_max[e][1]:
                    min_max[e][1] = max_value
                    flag_u = True
                    length_ = 0
                    for iv in vars__:
                        for v in iv:
                            length_ = length_ + v[1]
                    if length_ > length:
                        length = length_
                else:
                    max_value = min_max[e][1]
            
            query = []
            sum = 0
            for ii in e:
                sum = sum + model.get_initial_state()[ii]
            for ii in range(N + 1):
                for iii in vars__[ii]:
                    for iv in e:
                        sum = sum + (Int(iii[2]) * model.get_reactions_vector()[iii[0]][iv])
                query.append(sum > max_value)
            solver.add(Or(query))
        solver.pop()

            #
    #lower-bound
    flag_l = False
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars__ = [None] * (N + 1)
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars__[state_index] == None:
                        vars__[state_index] = []
                        vars__[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars__[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            min_value = 0
            for ii in e:
                min_value = min_value + model.get_initial_state()[ii]
            for ii in range(N + 1):
                for iii in vars__[ii]:
                    for iv in e:
                        min_value = min_value + (model.get_reactions_vector()[iii[0]][iv] * iii[1])
                if min_value < min_max[e][0]:
                    min_max[e][0] = min_value
                    flag_l = True
                    length_ = 0
                    for iv in vars__:
                        for v in iv:
                            length_ = length_ + v[1]
                    if length_ > length:
                        length = length_
                else:
                    min_value = min_max[e][0]
            
            query = []
            sum = 0
            for ii in e:
                sum = sum + model.get_initial_state()[ii]
            for ii in range(N + 1):
                for iii in vars__[ii]:
                    for iv in e:
                        sum = sum + (Int(iii[2]) * model.get_reactions_vector()[iii[0]][iv])
                query.append(sum < min_value)
            solver.add(Or(query))
        solver.pop()
    return (flag_u or flag_l, length, min_max)


################################################################################################
################################################################################################
################################################################################################

def avg_prob(model_name):
    # rate_constants = model.get_reaction_rate_constants()
    # sum = 0
    # for c in rate_constants:
    #     sum = sum + c
    # for i, e in enumerate(rate_constants):
    #     rate_constants[i] = rate_constants[i] / sum
    # return rate_constants
    #

    if "enzym" in model_name:
        dist =  [0.25095102429706134, 0.2262466319108741, 0.022400865227274323, 0.25135835531906425, 0.22659777934363523, 0.02244534390209073]
        policy = [1.0, 1.0, 0.5, 1.0, 1.0, 2.0]
        b = [None] * len(dist)
        for i, e in enumerate(b):
            b[i] = dist[i] * policy[i]
        return b
    elif "motil" in model_name:
        dist = [0.052953890489913544, 0.0011321531494442158, 0.3866817620419926, 0.000977768629065459, 0.14203375874845617, 0.0007719226018937835, 0.11645738987237546, 0.013894606834088103, 0.08100041169205434, 0.06329765335529024, 0.07889048991354466, 0.06190819267188143]
        policy = [3.33, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.3, 3.33, 0.3, 3.33]
        b = [None] * len(dist)
        for i, e in enumerate(b):
            b[i] = dist[i] * policy[i]
        return b
    elif "yeast" in model_name:
        dist = [0.0002426220054642695, 0.0006821546240589607, 0.14771460619635504, 0.011459494432000788, 0.2850773401595679, 0.16413906108800147, 0.16413906108800147, 0.22654566040655008]
        policy = [1.37, 0.725, 1.37, 0.725, 1.37, 0.725, 0.725, 1.37]
        b = [None] * len(dist)
        for i, e in enumerate(b):
            b[i] = dist[i] * policy[i]
        return dist
    elif "circuit" in model_name:
        dist =  [0.16974848319588945, 0.18720692110609946, 0.00026414135479406513, 0.05060064538872468, 0.00026213267528992775, 0.049973937383433815, 0.20472963676044187, 0.18632711948328728, 0.04995585926789658, 0.05030637384136855, 0.04977909547153249, 1.3056416776892954E-05, 0.00037964042628196436, 0.00023702418148821057, 0.00021593304669476808]
        policy = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.66, 1.66, 0.6]
        b = [None] * len(dist)
        for i, e in enumerate(b):
            b[i] = dist[i] * policy[i]
        return b
    else:
        raise Exception("Average probability vector for model is not defined")