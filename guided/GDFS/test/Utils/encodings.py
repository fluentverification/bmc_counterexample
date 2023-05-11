from z3 import Int, And, Or, simplify

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

# IMPORTANT: check if this encoding allows for 
# reactions of type s1 -> s1 + s2 to happen even
# when s1==0
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
            r_constraints.append(y >= 0)
        reaction_var_name = 'selected_reaction.' + str(bound-1) 
        selected_reaction = Int(reaction_var_name)
        r_constraints.append(selected_reaction == (j))
        constraints.append(And(r_constraints))
    return simplify(Or(constraints))

def target_encoding(bound, target_var, target_value):
    x = Int(target_var + '.' + str(bound))
    target_constraint = (x==target_value)
    return target_constraint