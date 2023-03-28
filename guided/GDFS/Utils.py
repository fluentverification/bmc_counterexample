from z3 import *
from Models import enzymatic_futile_cycle, motility_regulation, yeast_polarization
import itertools

def get_total_outgoing_rate(var_values, model):
    rate = 0
    for i, r in enumerate(model.get_reactions_vector()):
        comb = 1
        for j, coefficient in enumerate(r[0]):
            if coefficient>0:
                for c in range(coefficient):
                    comb = comb * var_values[j]
        rate = rate + model.get_reaction_rates()[i] * comb
    return rate

def get_reaction_rate(var_values, model, reaction):
    rate = model.get_reaction_rates()[reaction]
    for i, coefficient in enumerate(model.get_reactions_vector()[reaction][0]):
        if coefficient>0:
            for c in range(coefficient):
                rate = rate * var_values[i]
    return rate

def is_target(var_values, target_index, target_value):
    if var_values[target_index] == target_value:
        return True
    return False

#returns all the subsets of a tuple with cardinality in range [L, U]
#returns a list of tuples where each tuple is a subset
def get_subsets (tuple, L, U):
    if U>len(tuple):
        U = len(tuple)
    return list(itertools.chain.from_iterable(itertools.combinations(tuple, r) for r in range(L, U+1)))

def get_min_max(model, bound, target_index, target_value, maximum_comb):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in get_subsets(index_tuple, 1, maximum_comb):
        min_max[s] = [-1, -1]
    
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
    #third constraint (reaching the target)
    vars = []
    for i, r in enumerate(model.get_reactions_vector()):
        if r[2][target_index]!=0:
            vars.append([Int("n_i_" + str(i)), r[2][target_index]])
            vars.append([Int("n_t_" + str(i)), r[2][target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    #fourth constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[2][i]!=0:
                vars.append([Int("n_i_"+str(j)), r[2][i]])
                vars.append([Int("n_t_"+str(j)), r[2][i]])
                vars_i.append([Int("n_i_"+str(j)), r[2][i]])
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
                    max_value = max_value + (model.get_reactions_vector()[v[0]][2][j] * v[1])
            min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][2][j])
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
                    min_value = min_value + (model.get_reactions_vector()[v[0]][2][j] * v[1])
            min_max[e][0] = min_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][2][j])
            solver.add(sum<min_value)
        solver.pop()


    return min_max

def check_sat(var_values, min_max):
    for e in min_max:
        value = 0
        for j in e:
            value = value + var_values[j]
        if value > min_max[e][1] or value < min_max[e][0]:
            return False
    return True


target_index = 5
target_value = 50
#original bound is 100
bound = 103
model = yeast_polarization()
#print('here')
print(get_min_max(model = model, bound = bound, target_index=target_index, target_value=target_value, maximum_comb=100))
# # print(get_subsets((1,2,3,4), 1, 4))