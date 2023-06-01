from z3 import Int, And, Or, simplify, sat

def initial_state_encoding(model):
    species_vector = model.get_species_tuple()
    initial_state = model.get_initial_state()
    constraints = []
    for i, s in enumerate(species_vector): 
        var_name = s + '.' + '0'
        x = Int(var_name)
        init_value = initial_state[i]
        constraints.append(x == init_value)
    return And(constraints)

def state_encoding(model, state):
    species_vector = model.get_species_tuple()
    constraints = []
    for i, s in enumerate(species_vector): 
        var_name = s + '.' + '0'
        x = Int(var_name)
        var_value = state.var_values[i]
        constraints.append(x == var_value)
    return And(constraints)


# IMPORTANT: check if this encoding allows for 
# reactions of type s1 -> s1 + s2 to happen even
# when s1==0
# IMPORTANT: added new encoding by forcing reactions_lhs_vector.
# Check to see if get_reactions_lhs_vector works for all models.
def bound_encoding(model, bound):
    species_vector = model.get_species_tuple()
    constraints = []
    for j, r in enumerate(model.get_reactions_vector()): 
        r_constraints = []
        for i, s in enumerate(species_vector): 
            var_name_prev = s + '.' + str(bound-1)
            var_name_curr = s + '.' + str(bound)
            x = Int(var_name_prev)
            y = Int(var_name_curr)
            r_constraints.append(y == (x + r[i]))
            #r_constraints.append(y >= 0)
            r_constraints.append(x >= model.get_reactions_lhs_vector()[j][i])
        reaction_var_name = 'selected_reaction.' + str(bound-1) 
        selected_reaction = Int(reaction_var_name)
        r_constraints.append(selected_reaction == (j))
        constraints.append(And(r_constraints))
    return simplify(Or(constraints))

#
def target_encoding(bound, target_var, target_value):
    x = Int(target_var + '.' + str(bound))
    target_constraint = (x==target_value)
    return target_constraint
#

# target_encoding that accepts all the states in the graph as 
# well as the targets
def target_encoding_2(model, bound, target_var, target_value, diag):
    target_const_temp = []
    for n in diag.nodes.values():
        var_values = []
        for i, s in enumerate(model.get_species_tuple()):
            x = Int (s + '.' + str(bound))
            temp_const = (x == n.var_values[i])
            var_values.append(temp_const)
        target_const_temp.append(And(var_values))
    target_const_temp = Or(target_const_temp)
    x = Int(target_var + '.' + str(bound))
    property_constraint = (x==target_value)
    target_const = Or(target_const_temp, property_constraint)
    return target_const
#


# new unsat encoding.
def check_sat(model, solver, current_state, bound, target_index, target_value):
    solver.pop()
    solver.push()
    vars = []
    for i,r in enumerate (model.get_reactions_vector()):
        x = Int("n" + str(i))
        vars.append(x)
    
    constraints = []
    #first constraint (n1,n2,...>=0)
    for i in vars:
        constraints.append(i>=0)
    #

    #second constraint (n1+n2+...=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum==bound))
    #
    
    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            vars.append([Int("n" + str(i)), r[target_index]])
    sum = current_state[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    #
    
    #fourth constraint (species population >=0)
    for i, s in enumerate(model.get_species_tuple()):
        if i==target_index:
            continue
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n" + str(j)), r[i]])
        sum = current_state[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0) 

    solver.add(And(constraints))
    if (solver.check()==sat):
        return True
    else:
        return False


# new unsat encoding. Accepts all the nodes in the diag
def check_sat_2(model, diag, solver, current_state, bound, target_index, target_value):
    return True
    quit()
    solver.pop()
    solver.push()
    vars = []
    for i,r in enumerate (model.get_reactions_vector()):
        x = Int("n" + str(i))
        vars.append(x)
    
    constraints = []
    #first constraint (n1,n2,...>=0)
    for i in vars:
        constraints.append(i>=0)
    #

    #second constraint (n1+n2+...=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum==bound))
    #
 
    #third constraint (reaching the target)
    target_const = []
    for n in diag.nodes.values():
        target_const_temp = []
        for i, s in enumerate(model.get_species_tuple()):
            vars = []
            for j, r in enumerate(model.get_reactions_vector()):
                vars.append([Int("n" + str(j)), r[i]])
            sum = current_state[j]
            for j in vars:
                sum = sum + j[0]*j[1]
            target_const_temp.append(sum== n.var_values[i])
        target_const.append(And(target_const_temp))
    target_const = Or(target_const)

    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[target_index]!=0:
            vars.append([Int("n" + str(i)), r[target_index]])
    sum = current_state[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]

    constraints.append(Or((sum==target_value), target_const))
    #
    
    #fourth constraint (species population >=0)
    for i, s in enumerate(model.get_species_tuple()):
        if i==target_index:
            continue
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n" + str(j)), r[i]])
        sum = current_state[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0) 

    solver.add(And(constraints))
    if (solver.check()==sat):
        return True
    else:
        return False