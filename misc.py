def get_min_max_bound_only(model, bound, target_index, target_value, subsets, min_max_prev):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_i_" + str(i))
        vars.append(x)
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_t_" + str(i))
        vars.append(x)
    constraints = []
    
    #first constraint (n_i_0,n_t_0,n_i_1,n_t_1,... >=0)
    for i in vars:
        constraints.append(i>=0)
    
    #second constraint (n_i_0+n_i_t+n_i_1+...<=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum<=bound))

    #third constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n_i_"+str(j)), r[i]])
                vars.append([Int("n_t_"+str(j)), r[i]])
                vars_i.append([Int("n_i_"+str(j)), r[i]])
        sum = model.get_initial_state()[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)
        sum = model.get_initial_state()[i]
        for j in vars_i:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)

    solver = Solver()
    solver.add(And(constraints))

    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            max_value = 0
            for j in e:
                max_value = max_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    max_value = max_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            if max_value > min_max[e][1]:
                min_max[e][1] = max_value
            elif min_max[e][1] != -1:
                max_value = min_max[e][1]
            else:
                min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum>max_value)
        solver.pop()
    
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            min_value = 0
            for j in e:
                min_value = min_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    min_value = min_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][0] = min_value
            if min_value < min_max[e][0]:
                min_max[e][0] = min_value
            elif min_max[e][0] != -1:
                min_value = min_max[e][0]
            else:
                min_max[e][0] = min_value

            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum<min_value)
        solver.pop()
    return min_max
#

def get_min_max_probability_jan5th(model, prob_thresh, target_index, target_value, subsets, min_max_prev):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_i_" + str(i))
        vars.append(x)
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_t_" + str(i))
        vars.append(x)
    constraints = []
    
    #first constraint (n_i_0,n_t_0,n_i_1,n_t_1,... >=0)
    for i in vars:
        constraints.append(i>=0)
    
    #second constraint (probability threshold)
    avg_population = [None] * len(model.get_species_tuple())
    for i, s in enumerate(model.get_species_tuple()):
        pop = model.get_initial_state()[i]
        change = 0
        for j in range (0, len(vars)):
            change = change + (vars[j] * model.get_reactions_vector()[j%(len(model.get_reactions_vector()))][i])
        sum_vars = 0
        for k in vars:
            sum_vars = sum_vars + k
        pop = pop + (change / (sum_vars + 1))
        avg_population[i] = pop

    avg_rate = [None] * len(model.get_reactions_vector())
    for i, r in enumerate(model.get_reactions_vector()):
        avg_rate[i] = model.get_reaction_rate_constants()[i]
        for j, s in enumerate(model.reactions_lhs_vector[i]):
            if s > 0:
                for k in range (s):
                    avg_rate[i] = avg_rate[i] * avg_population[j]
    total_rate = 0
    for v in avg_rate:
        total_rate = total_rate + v
    avg_prob = [None] * len(model.get_reactions_vector())
    for i, v in enumerate(avg_rate):
        avg_prob[i] = v/total_rate
    
    prob = 0
    for i, _ in enumerate(model.get_reactions_vector()):
        # prob = prob * (avg_prob[i] ** (vars[i] + vars[len(model.get_reactions_vector()) + i]))
        # prob = prob * avg_prob[i]
        prob = prob + (((avg_prob[i]-1) - 0.5*((avg_prob[i]-1)**2))  * (vars[i] + vars[len(model.get_reactions_vector()) + i]))

    constraints.append((prob>prob_thresh))
    

    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            vars.append([Int("n_i_" + str(i)), r[target_index]])
            vars.append([Int("n_t_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    
    #fourth constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n_i_"+str(j)), r[i]])
                vars.append([Int("n_t_"+str(j)), r[i]])
                vars_i.append([Int("n_i_"+str(j)), r[i]])
        sum = model.get_initial_state()[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)
        sum = model.get_initial_state()[i]
        for j in vars_i:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)

    
    
    #
    #
    solver = Solver()
    solver.add(And(constraints))


    flag1 = False
    #print('hereee')
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            print('here')
            flag1 = True
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            max_value = 0
            for j in e:
                max_value = max_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    max_value = max_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            if max_value > min_max[e][1]:
                min_max[e][1] = max_value
            elif min_max[e][1] != -1:
                max_value = min_max[e][1]
            else:
                min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum>max_value)
        solver.pop()
    
    flag2 = False
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            flag2 = True
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            min_value = 0
            for j in e:
                min_value = min_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    min_value = min_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][0] = min_value
            if min_value < min_max[e][0]:
                min_max[e][0] = min_value
            elif min_max[e][0] != -1:
                min_value = min_max[e][0]
            else:
                min_max[e][0] = min_value

            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum<min_value)
        solver.pop()
    return (flag1 and flag2, min_max)
#

def get_min_max_probability(model, prob_thresh, target_index, target_value, subsets, min_max_prev):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in subsets:
        min_max[s] = min_max_prev[s]
    
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_i_" + str(i))
        vars.append(x)
    for i, r in enumerate(model.get_reactions_vector()):
        x= Int("n_t_" + str(i))
        vars.append(x)
    constraints = []
    
    #first constraint (n_i_0,n_t_0,n_i_1,n_t_1,... >=0)
    for i in vars:
        constraints.append(i>=0)
    
    #second constraint (probability threshold)
    avg_population = [None] * len(model.get_species_tuple())
    for i, s in enumerate(model.get_species_tuple()):
        pop = model.get_initial_state()[i]
        change = 0
        for j in range (0, len(vars)):
            change = change + (vars[j] * model.get_reactions_vector()[j%(len(model.get_reactions_vector()))][i])
        sum_vars = 0
        for k in vars:
            sum_vars = sum_vars + k
        pop = pop + (change / (sum_vars + 1))
        avg_population[i] = pop

    avg_rate = [None] * len(model.get_reactions_vector())
    for i, r in enumerate(model.get_reactions_vector()):
        avg_rate[i] = model.get_reaction_rate_constants()[i]
        for j, s in enumerate(model.reactions_lhs_vector[i]):
            if s > 0:
                for k in range (s):
                    avg_rate[i] = avg_rate[i] * avg_population[j]
    total_rate = 0
    for v in avg_rate:
        total_rate = total_rate + v
    avg_prob = [None] * len(model.get_reactions_vector())
    for i, v in enumerate(avg_rate):
        avg_prob[i] = v/total_rate
    
    prob = 0
    for i, _ in enumerate(model.get_reactions_vector()):
        prob = prob + (math.log(model.get_reaction_rate_constants()[i]) * (vars[i] + vars[len(model.get_reactions_vector()) + i]))
        # prob = prob + (avg_prob[i] / len(model.get_reactions_vector()))
    constraints.append((prob>prob_thresh))
    

    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            vars.append([Int("n_i_" + str(i)), r[target_index]])
            vars.append([Int("n_t_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    

    #fourth constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n_i_"+str(j)), r[i]])
                vars.append([Int("n_t_"+str(j)), r[i]])
                vars_i.append([Int("n_i_"+str(j)), r[i]])
        sum = model.get_initial_state()[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)
        sum = model.get_initial_state()[i]
        for j in vars_i:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0)

    
    
    #
    #
    solver = Solver()
    solver.add(And(constraints))

    flag1 = False
    #print('hereee')
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            # print('here')
            flag1 = True
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            max_value = 0
            for j in e:
                max_value = max_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    max_value = max_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            if max_value > min_max[e][1]:
                min_max[e][1] = max_value
            elif min_max[e][1] != -1:
                max_value = min_max[e][1]
            else:
                min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum>max_value)
        solver.pop()
    
    flag2 = False
    for i, e in enumerate(min_max):
        solver.push()
        while (solver.check()==sat):
            flag2 = True
            assignment = solver.model()
            vars = []
            for d in assignment.decls():
                if "n_i_" in str(d.name()):
                    r_index = int(d.name()[4:])
                    vars.append([r_index, assignment[d].as_long(), d.name()])
            min_value = 0
            for j in e:
                min_value = min_value + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    min_value = min_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][0] = min_value
            if min_value < min_max[e][0]:
                min_max[e][0] = min_value
            elif min_max[e][0] != -1:
                min_value = min_max[e][0]
            else:
                min_max[e][0] = min_value

            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum<min_value)
        solver.pop()
    return (flag1 and flag2, min_max)
#







################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
######################## test test test test test test test test test te########################
def get_min_max_abs(model, prob_thresh, K, target_index, target_value, subsets, min_max_prev, N):
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
    species_pops = [None] * (N + 1)
    initial_valuation = list(model.get_initial_state())
    species_pops[0] = initial_valuation
    for i in range(1, N + 1):
        valuation = [None] * len(model.get_species_tuple())
        for ii, s in enumerate(model.get_species_tuple()):
            valuation[ii] = initial_valuation[ii]
            for iii, r in enumerate(model.get_reactions_vector()):
                x = Int("n_" + str(i-1) + "_" + str(iii))
                if r[ii] > 0:
                    valuation[ii] = valuation[ii] + (x * r[ii])
        species_pops[i] = valuation
        initial_valuation = valuation

    prob = 0
    for i in range(N+1):
        rate_vector = [None] * len(model.get_reactions_vector())
        for ii, r in enumerate(model.get_reactions_vector()):
            rate_vector[ii] = model.get_reaction_rate_constants()[ii]
            for iii, s in enumerate(model.reactions_lhs_vector[ii]):
                if s > 0:
                    for iv in range (s):
                        rate_vector[ii] = rate_vector[ii] * species_pops[i][iii]

        total_rate = 0
        for ii in rate_vector:
            total_rate = total_rate + ii
        
        prob_vector = [None] * len(model.get_reactions_vector())
        for ii, v in enumerate(rate_vector):
            prob_vector[ii] = v/total_rate
        
        for ii, _ in enumerate(model.get_reactions_vector()):
            x = Int("n_" + str(i) + "_" + str(ii))
            # prob = prob + (((prob_vector[ii]-1) - 0.5*((prob_vector[ii]-1)**2))  * x)
            prob = prob + (((prob_vector[ii]-1))  * x)


    #constraints.append(prob > prob_thresh)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append(sum <= K)
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

