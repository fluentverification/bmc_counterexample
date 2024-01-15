from z3 import Solver, Int, And, Or, sat, Real
import math
from Poisson import poisson_cdf, poisson_cdf_step, avg_discrete_probability


def get_min_max_reaction(model, model_name, prob_thresh, poisson_step, division_factor, subsets, min_max_prev):
    min_max = {}
    for s in subsets:
        min_max[s] = min_max_prev[s]

    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        x = Int("n_" + str(i))
        vars.append(x)

    constraints = []
    
    #first constraint: n_0, n_1, ... >=0
    for i in vars:
        constraints.append(i>=0)
    #
    
    #second constraint: probability threshold
    prob_vars = []
    if poisson_step == 1:
        prob_vector = poisson_cdf(model_name= model_name, cut_off=prob_thresh, division_factor=division_factor)
        if prob_vector == None:
            return True, None
    else:
        prob_vector = poisson_cdf_step(model_name=model_name, step=poisson_step, cut_off=prob_thresh, division_factor=division_factor)
    for i, _ in enumerate(model.get_reactions_vector()):
        x = Int("n_" + str(i))
        prob_dict = prob_vector[i]
        if prob_dict == None:
            constraints.append(x == 0)
        else:
            # constraints.append(x < list(prob_dict.keys())[-1][1])
            prob_approx = []
            prob_value = Real("prob_" + str(i))
            prob_vars.append(prob_value)
            for key, value in prob_dict.items():
                if poisson_step == 1:
                    prob_approx.append(And(x == key, prob_value == value))
                else:
                    prob_approx.append(And(x>=key[0], x<key[1], prob_value == value))
                    if key == list(prob_dict.keys())[-1]:
                        prob_approx.append(And(x==key[1], prob_value == value))
            constraints.append(Or(prob_approx))
    
    sum = 0
    for i in prob_vars:
        sum = sum + i
    constraints.append(sum > prob_thresh)
        
    #querying the solver
    solver = Solver()
    solver.add(And(constraints))
    flag_u = False
    for i, e in enumerate(min_max):
        solver.push()
        while(solver.check()==sat):
            assignment = solver.model()
            vars_ = []
            for d in assignment.decls():
                if "n_" in d.name():
                    r_index = int(d.name()[2:])
                    vars_.append([r_index, assignment[d].as_long(), d.name()])
            
            max_value = 0
            for i in e:
                for v in vars_:
                    if v[0] == i:
                        max_value = max_value + v[1]
            
            if max_value > min_max[e]:
                min_max[e] = max_value
                flag_u = True
            else:
                max_value = min_max[e]
            
            sum = 0
            for i in e:
                for v in vars_:
                    if v[0] == i:
                        sum = sum + Int(v[2])
            solver.add(sum > max_value)
        solver.pop()
    return(flag_u, min_max)

def get_min_max_species(model, min_max_reaction, subsets, min_max_prev, target_index, target_value):
    min_max = {}
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    steps = 0
    for key, value in min_max_reaction.items():
        if len(key) == 1:
            steps = steps + value
    
    steps_log = math.floor(math.log(steps))
    max_len = math.ceil(steps / steps_log)
    steps = steps_log

    vars = []
    for i in range(steps):
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            vars.append(x)
        
    constraints = []
    
    #first constraint : n_0_0,n_0_1,n_1_0,n_1_1,... >=0
    for i in vars:
        constraints.append(i>=0)
    #
        
    #second constraint : some of reactions in each step must be less than maximum length for each step
    for i in range(steps):
        sum = 0
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            sum = sum + x
        constraints.append(sum<=max_len)
    #
    
    #third constraint : species population must be greater than zero at each arbitrary state
    species_pops = [None] * steps
    initial_valuation = list(model.get_initial_state())
    for i in range(0, steps):
        valuation = [None] * len(model.get_species_tuple())
        for s, _ in enumerate(model.get_species_tuple()):
            valuation[s] = initial_valuation[s]
            for j, r in enumerate(model.get_reactions_vector()):
                x = Int("n_" + str(i) + "_" + str(j))
                if r[s] != 0:
                    valuation[s] = valuation[s] + (x * r[s])
        species_pops[i] = valuation
        initial_valuation = valuation
        for s in valuation:
            constraints.append(s>=0)
    #
            
    #fourth constraint : conforming to constraints on reactions
    for key, value in min_max_reaction.items():
        sum = 0
        for r in key:
            for i in range(steps):
                x = Int("n_" + str(i) + "_" + str(r))
                sum = sum + x
        constraints.append(sum <= value)
    #
    
    #fifth constraint : reaching the target
    vars_ = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            for ii in range(steps):
                vars_.append([Int("n_" + str(ii) + "_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars_:
        sum = sum + i[0]*i[1]

    constraints.append(sum==target_value)   
    #
        
        
    #querying the solver
    solver = Solver()
    solver.add(And(constraints))

    #upper-bound
    flag_u = False
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars_ = [None] * steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_[state_index] == None:
                        vars_[state_index] = []
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            pop = 0
            max_value = 0
            for s in key:
                pop = pop + model.get_initial_state()[s]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        pop = pop + (model.get_reactions_vector()[v[0]][s] * v[1])
                if pop > min_max[key][1]:
                    min_max[key][1] = pop
                    max_value = pop
                    flag_u = True
                else:
                    max_value = min_max[key][1]
            
            query = []
            sum = 0
            for i in key:
                sum = sum + model.get_initial_state()[i]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum > max_value)
            solver.add(Or(query))
        solver.pop()

            #
    #lower-bound
    flag_l = False
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars_ = [None] * steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_[state_index] == None:
                        vars_[state_index] = []
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            pop = 0
            min_value = 0
            for s in key:
                pop = pop + model.get_initial_state()[s]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        pop = pop + (model.get_reactions_vector()[v[0]][s] * v[1])
                if pop < min_max[key][0]:
                    min_max[key][0] = pop
                    min_value = pop
                    flag_l = True
                else:
                    min_value = min_max[key][0]
            
            query = []
            sum = 0
            for i in key:
                sum = sum + model.get_initial_state()[i]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum < min_value)
            solver.add(Or(query))
        solver.pop()
    
    return (flag_u or flag_l, min_max)

def get_min_max(model, model_name, prob_thresh, max_len, division_factor,subsets, min_max_prev, target_index, target_value):
    min_max = {}
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    steps = math.floor(math.log(max_len, 10))
    max_segment_len = math.ceil(max_len / steps)

    vars = []
    for i in range(steps):
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            vars.append(x)
        
    constraints = []
    
    #first constraint : n_0_0,n_0_1,n_1_0,n_1_1,... >=0
    for i in vars:
        constraints.append(i>=0)
    #
        
    #second constraint : some of reactions in each step must be less than maximum length for each step
    for i in range(steps):
        sum = 0
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            sum = sum + x
        constraints.append(sum<=max_len)
    #
        
    #third constraint : species population must be greater than zero at each arbitrary state
    species_pops = [None] * steps
    initial_valuation = list(model.get_initial_state())
    for i in range(0, steps):
        valuation = [None] * len(model.get_species_tuple())
        for s, _ in enumerate(model.get_species_tuple()):
            valuation[s] = initial_valuation[s]
            for j, r in enumerate(model.get_reactions_vector()):
                x = Int("n_" + str(i) + "_" + str(j))
                if r[s] != 0:
                    valuation[s] = valuation[s] + (x * r[s])
        species_pops[i] = valuation
        initial_valuation = valuation
        for s in valuation:
            constraints.append(s>=0)
    #
    
    #fourth constraint (probability threshold)
    prob = 0
    #reachability_probability
    avg_prob_ = avg_discrete_probability(model_name)
    for i, _ in enumerate(model.get_reactions_vector()):
        sum = 0
        for ii in range(steps):
            x = Int("n_" + str(ii) + "_" + str(i))
            sum = sum + x
        prob = prob + (math.log(avg_prob_[i], division_factor) * sum)

    constraints.append((prob>prob_thresh))
    #probability of those reactions firing in allotted time
    # prob_vars = []
    # prob_vector = poisson_cdf(model_name= model_name, cut_off=prob_thresh, division_factor=division_factor)
    # for i, _ in enumerate(model.get_reactions_vector()):
    #     x = 0
    #     for j in range(steps):
    #         x_ = Int("n_" + str(j) + "_" + str(i))
    #         x = x + x_
    #     prob_dict = prob_vector[i]
    #     if prob_dict == None:
    #         constraints.append(x == 0)
    #     else:
    #         # constraints.append(x < list(prob_dict.keys())[-1][1])
    #         prob_approx = []
    #         prob_value = Real("prob_" + str(i))
    #         prob_vars.append(prob_value)
    #         for key, value in prob_dict.items():
    #             prob_approx.append(And(x == key, prob_value == value))
    #         constraints.append(Or(prob_approx))
    
    # sum = 0
    # for i in prob_vars:
    #     sum = sum + i
    # constraints.append(sum + prob > prob_thresh)
    #

    #fifth constraint : reaching the target
    vars_ = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            for ii in range(steps):
                vars_.append([Int("n_" + str(ii) + "_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars_:
        sum = sum + i[0]*i[1]

    # constraints.append(sum==target_value) 


    #querying the solver
    solver = Solver()
    solver.add(And(constraints))

    #upper-bound
    flag_u = False
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars_ = [None] * steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_[state_index] == None:
                        vars_[state_index] = []
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            pop = 0
            max_value = 0
            for s in key:
                pop = pop + model.get_initial_state()[s]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        pop = pop + (model.get_reactions_vector()[v[0]][s] * v[1])
                if pop > min_max[key][1]:
                    min_max[key][1] = pop
                    max_value = pop
                    flag_u = True
                else:
                    max_value = min_max[key][1]
            
            query = []
            sum = 0
            for i in key:
                sum = sum + model.get_initial_state()[i]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum > max_value)
            solver.add(Or(query))
        solver.pop()

            #
    #lower-bound
    flag_l = False
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars_ = [None] * steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    state_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_[state_index] == None:
                        vars_[state_index] = []
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_[state_index].append([r_index, assignment[d].as_long(), d.name()])
            
            pop = 0
            min_value = 0
            for s in key:
                pop = pop + model.get_initial_state()[s]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        pop = pop + (model.get_reactions_vector()[v[0]][s] * v[1])
                if pop < min_max[key][0]:
                    min_max[key][0] = pop
                    min_value = pop
                    flag_l = True
                else:
                    min_value = min_max[key][0]
            
            query = []
            sum = 0
            for i in key:
                sum = sum + model.get_initial_state()[i]
            for i in range(steps):
                for v in vars_[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum < min_value)
            solver.add(Or(query))
        solver.pop()
    
    return (flag_u or flag_l, min_max)