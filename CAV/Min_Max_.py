import math
from z3 import Int, Real, And, Or, Solver, sat, simplify
from Probability import t_abstract_prob, poisson_at_least_k

def get_min_max_species(model, 
                        model_name, 
                        prob_thresh, 
                        num_steps, 
                        division_factor,
                        poisson_step,
                        subsets,
                        min_max_prev, 
                        target_index, 
                        target_value, 
                        lower_bound):
    #for 3 steps we have
    #initial_state (--step_0-->) arbitrary_state_0 (--step_1-->) arbitrary_state_1 (--step_2-->) target_state
    #

    #constraints fed into the solver
    constraints = []
    #

    #initializing the min_max 
    min_max = {}
    for k, e in min_max_prev.items():
        min_max[k] = e
    #
    
    #declaring the variables
    vars = []
    for i in range(num_steps):
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            vars.append(x)
    #
    
    #the number of steps (arbitrary states + 1) the trace takes to reach the target and the 
    #maximum number of reactions that could be taken within each step
    sum = 0
    for i in vars:
        sum = sum + i
    max_segment_len = (sum / num_steps)
    #
    
    #first constraint : n_0_0,n_0_1,n_1_0,n_1_1,... >=0
    for i in vars:
        constraints.append(i>=0)
    #
    
    #second constraint : sum of the reactions in each step must be less than maximum length for each step
    for i in range(num_steps):
        sum = 0
        for ii, _ in enumerate(model.get_reactions_vector()):
            x= Int("n_" + str(i) + "_" + str(ii))
            sum = sum + x
        constraints.append(sum<=(max_segment_len + 1))
    #
    
    #third constraint : there must be enough reactant molecules for the reaction firings in each step
    initial_population = list(model.get_initial_state())
    for i in range(num_steps):
        necessary_reactants = [0] * len(model.get_species_tuple())
        for ii, lhs in enumerate(model.reactions_lhs_vector):
            x = Int("n_" + str(i) + "_" + str(ii))
            for s, _ in enumerate(model.get_species_tuple()):
                necessary_reactants[s] = necessary_reactants[s] + (x * lhs[s])

        available_reactants = [0] * len(model.get_species_tuple())
        for ii, s in enumerate(initial_population):
            available_reactants[ii] = s
        for ii, r in enumerate(model.get_reactions_vector()):
            x = Int("n_" + str(i) + "_" + str(ii))
            for s, _ in enumerate(model.get_species_tuple()):
                if r[s] > 0:
                    available_reactants[s] = available_reactants[s] + (x * r[s])
        
        for ii, s in enumerate(necessary_reactants):
            constraints.append(s<=available_reactants[ii])

        for ii, s in enumerate(initial_population):
            sum = s
            for iii, r in enumerate(model.get_reactions_vector()):
                x = Int("n_" + str(i) + "_" + str(iii))
                sum = sum + (r[ii] * x)
            initial_population[ii] = sum
    #
    
    #fourth constraint : estimated probability of the path must be greater than probability threshold
    trace_prob = 0
    t_abs_probs = t_abstract_prob(model_name=model_name, division_factor=division_factor)
    #time abstract probability
    for i in range(len(model.get_reactions_vector())):
        sum = 0
        for ii in range(num_steps):
            x = Int("n_" + str(ii) + "_" + str(i))
            sum = sum + x
        trace_prob = trace_prob + (sum * t_abs_probs[i])
    constraints.append(trace_prob >= prob_thresh)
    #

    #fifth constraint : reaching the target
    if lower_bound:
        vars_ = []
        for i, r in enumerate(model.get_reactions_vector()):
            if r[target_index]!=0:
                for ii in range(num_steps):
                    vars_.append([Int("n_" + str(ii) + "_" + str(i)), r[target_index]])
        sum = model.get_initial_state()[target_index]
        for i in vars_:
            sum = sum + i[0]*i[1]
        constraints.append(sum==target_value) 
    #
        
    #free variables representing an arbitrary state among each trace segment (step)
    for i in range(num_steps):
        for ii in range(len(model.get_reactions_vector())):
            x = Int("i_" + str(i) + "_" + str(ii))
            y = Int("n_" + str(i) + "_" + str(ii))
            constraints.append(And(x>=0, x<=y))
    #
    
    #Calling the solver to get lower and upper bounds
    return get_bounds_(constraints, min_max, model, num_steps)
    #

def get_bounds_(constraints, min_max, model, num_steps):
    solver = Solver()
    solver.add(And(constraints))
    flag = False
    max_len = 0

    #upper-bound
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            #vars stores n_0_0, n_0_1, n_1_0, n_1_1,...
            #vars_i stores i_0_0, i_i_1, i_1_0, i_1_1,...
            vars = [None] * num_steps
            vars_i = [None] * num_steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars[step_index] == None:
                        vars[step_index] = []
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                if "i_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_i[step_index] == None:
                        vars_i[step_index] = []
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])

            #what's the length of the modeled trace
            length = 0
            for i in range(num_steps):
                for v in vars[i]:
                    length = length + v[1]
            if length > max_len:
                max_len = length
                    
            #what would be the maximum value of the key for the current assignment
            max_value = 0
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                population_temp = initial_population
                for v in vars_i[i]:
                    for s in key:
                        population_temp = population_temp + (model.get_reactions_vector()[v[0]][s] * v[1])
                if population_temp > min_max[key][1]:
                    flag = True
                    max_value = population_temp
                    min_max[key][1] = population_temp
                else:
                    max_value = min_max[key][1]
                for s in key:
                    for v in vars[i]:
                        initial_population = initial_population + (model.get_reactions_vector()[v[0]][s] * v[1])
            
            #querying the solver to return a larger max value
            query = []
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                sum = initial_population
                for v in vars_i[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum > max_value)
                for v in vars[i]:
                    for s in key:
                        initial_population = initial_population + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])

            solver.add(Or(query))
        solver.pop()

    #lower-bound
    for key, value in min_max.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            #vars stores n_0_0, n_0_1, n_1_0, n_1_1,...
            #vars_i stores i_0_0, i_i_1, i_1_0, i_1_1,...
            vars = [None] * num_steps
            vars_i = [None] * num_steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars[step_index] == None:
                        vars[step_index] = []
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                if "i_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_i[step_index] == None:
                        vars_i[step_index] = []
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])

            #what's the length of the modeled trace
            length = 0
            for i in range(num_steps):
                for v in vars[i]:
                    length = length + v[1]
            if length > max_len:
                max_len = length
                    
            #what would be the minimum value of the key for the current assignment
            min_value = 0
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                population_temp = initial_population
                for v in vars_i[i]:
                    for s in key:
                        population_temp = population_temp + (model.get_reactions_vector()[v[0]][s] * v[1])
                if (population_temp < min_max[key][0]) and (population_temp >= 0):
                    flag = True
                    min_value = population_temp
                    min_max[key][0] = population_temp
                else:
                    min_value = min_max[key][0]
                for s in key:
                    for v in vars[i]:
                        initial_population = initial_population + (model.get_reactions_vector()[v[0]][s] * v[1])
            
            #querying the solver to return a smaller max value
            query = []
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                sum = initial_population
                for v in vars_i[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(And(sum < min_value, sum >= 0))
                for v in vars[i]:
                    for s in key:
                        initial_population = initial_population + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])

            solver.add(Or(query))
        solver.pop()
    return flag, max_len, min_max

def get_bounds(constraints, min_max, model, num_steps):
    solver = Solver()
    solver.add(And(constraints))
    flag = False
    max_len = 0

    #upper-bound
    items = list(min_max.values())[0]
    for key, _ in items.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            #vars stores n_0_0, n_0_1, n_1_0, n_1_1,...
            #vars_i stores i_0_0, i_i_1, i_1_0, i_1_1,...
            vars = [None] * num_steps
            vars_i = [None] * num_steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars[step_index] == None:
                        vars[step_index] = []
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                if "i_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_i[step_index] == None:
                        vars_i[step_index] = []
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])

            #what's the length of the modeled trace
            length = 0
            for i in range(num_steps):
                for v in vars[i]:
                    length = length + v[1]
            if length > max_len:
                max_len = length
                    
            #what would be the maximum value of the key for the current assignment
            max_value = [0] * num_steps
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                population_temp = initial_population
                for v in vars_i[i]:
                    for s in key:
                        population_temp = population_temp + (model.get_reactions_vector()[v[0]][s] * v[1])
                if population_temp > min_max[i][key][1]:
                    flag = True
                    max_value[i] = population_temp
                    min_max[i][key][1] = population_temp
                else:
                    max_value[i] = min_max[i][key][1]
                for s in key:
                    for v in vars[i]:
                        initial_population = initial_population + (model.get_reactions_vector()[v[0]][s] * v[1])
            
            #querying the solver to return a larger max value
            query = []
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                sum = initial_population
                for v in vars_i[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(sum > max_value[i])
                for v in vars[i]:
                    for s in key:
                        initial_population = initial_population + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
            solver.add(Or(query))
        solver.pop()

    for key, _ in items.items():
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            #vars stores n_0_0, n_0_1, n_1_0, n_1_1,...
            #vars_i stores i_0_0, i_i_1, i_1_0, i_1_1,...
            vars = [None] * num_steps
            vars_i = [None] * num_steps
            for d in assignment.decls():
                if "n_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars[step_index] == None:
                        vars[step_index] = []
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars[step_index].append([r_index, assignment[d].as_long(), d.name()])
                if "i_" in d.name():
                    first_occurrence = d.name().find("_")
                    second_occurrence = d.name().find("_", first_occurrence + 1)
                    step_index = int(d.name()[first_occurrence + 1:second_occurrence])
                    r_index = int(d.name()[second_occurrence + 1:])
                    if vars_i[step_index] == None:
                        vars_i[step_index] = []
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])
                    else:
                        vars_i[step_index].append([r_index, assignment[d].as_long(), d.name()])

            #what's the length of the modeled trace
            length = 0
            for i in range(num_steps):
                for v in vars[i]:
                    length = length + v[1]
            if length > max_len:
                max_len = length
                    
            #what would be the maximum value of the key for the current assignment
            min_value = [0] * num_steps
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                population_temp = initial_population
                for v in vars_i[i]:
                    for s in key:
                        population_temp = population_temp + (model.get_reactions_vector()[v[0]][s] * v[1])
                if (population_temp < min_max[i][key][0]) and (population_temp>=0):
                    flag = True
                    min_value[i] = population_temp
                    min_max[i][key][0] = population_temp
                else:
                    min_value[i] = min_max[i][key][0]
                for s in key:
                    for v in vars[i]:
                        initial_population = initial_population + (model.get_reactions_vector()[v[0]][s] * v[1])
            
            #querying the solver to return a larger max value
            query = []
            initial_population = 0
            for s in key:
                initial_population = initial_population + model.get_initial_state()[s]
            for i in range(num_steps):
                sum = initial_population
                for v in vars_i[i]:
                    for s in key:
                        sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
                query.append(And(sum < min_value[i], sum>=0))
                for v in vars[i]:
                    for s in key:
                        initial_population = initial_population + (Int(v[2]) * model.get_reactions_vector()[v[0]][s])
            solver.add(Or(query))
        solver.pop()
    return flag, max_len, min_max