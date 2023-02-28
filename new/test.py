from z3 import *
import Models
import numpy as np

def check_sat(model, state_vector, bound, target_var, target_value):
    vars = []
    for r in model.reactions_dict():
        x = Int("n" + str(r))
        vars.append(x)
    
    constraints = []
    #first constraint (n1,n2,...>=0)
    for i in vars:
        constraints.append(i>=0)

    #second constraint (n1+n2+...=bound)
    sum = 0
    for i in vars:
        sum = sum + i
    constraints.append((sum==bound))
    
    #third constraint (reaching the target)
    vars = []
    for i, s in enumerate(model.species_vector()):
        if s==target_var:
            target_index = i
    for i, r in enumerate(model.reactions_dict()):
        if model.reactions_dict()[r][2][target_index]!=0:
            vars.append([Int("n" + str(r)), model.reactions_dict()[r][2][target_index]])
    sum = state_vector[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    #fourth constraint (species population >=0)
    for i, s in enumerate(model.species_vector()):
        if s==target_var:
            continue
        vars = []
        for j, r in enumerate(model.reactions_dict()):
            if model.reactions_dict()[r][2][i]!=0:
                vars.append([Int("n" + str(r)), model.reactions_dict()[r][2][i]])
        sum = state_vector[i]
        for j in vars:
            sum = sum + j[0]*j[1]
        constraints.append(sum>=0) 

    solver = Solver()
    solver.add(And(constraints))
    if (solver.check()==sat):
        return True
    else:
        return False
           


constructor = getattr(Models, "yeast_polarization")
model = constructor()

state_vector = model.initial_state()

for i in range(120):
 print (str(i) + str(check_sat(model, state_vector, bound=i, target_var="Gbg", target_value=50)))

# def path_probability_encoding(L, U, P):
#     # define state variables S_i^j
#     S = [ [Bool(f"S_{i}^{j}") for j in range(num_states)] for i in range(U) ]
    
#     # define target state
#     theta = [Bool(f"theta_{j}") for j in range(num_states)]
    
#     # initial state encoding
#     enc_0 = get_encoding(0)
    
#     # path encoding up to bound L
#     path_encoding = And([get_encoding(i) for i in range(1, L+1)])
    
#     # target state is reached in at least L and at most U steps
#     target_reached = Or([S[i][-1] == theta[-1] for i in range(L, U+1)])
    
#     # for all states between L and U-1
#     # if the state was a target state, stop. All the following states would be the same.
#     # if not, ensure that the model takes an allowable transition in the next step.
#     between_L_U = And([Implies(S[i][-1] == theta[-1], And([S[j][-1] == S[i][-1] for j in range(i+1, U)])) 
#                       for i in range(L, U-1)])
#     between_L_U = And(between_L_U, 
#                       And([Implies(S[i][-1] != theta[-1], Implies(get_encoding(i+1), 
#                                                                  And([S[i+1][j] == S[i][j] for j in range(num_states)])))) 
#                       for i in range(L, U-1)])
    
#     # probability encoding
#     prob_encoding = Real(f"prob_encoding")
#     prob_path_reached = Real(f"prob_path_reached")
    
#     # set the probability of the initial state encoding
#     prob_encoding_init = If(enc_0, 1, 0)
    
#     # set the probability of the path encoding up to bound L
#     prob_path_encoding = 1
#     for i in range(1, L+1):
#         prob_path_encoding *= If(get_encoding(i), P, 1-P)
    
#     # set the probability that the target state is reached within L to U steps
#     prob_target_reached = 0
#     for i in range(L, U+1):
#         prob_target_reached += If(S[i][-1] == theta[-1], prob_path_encoding * (1-P)**(i-L) * P, 0)
    
#     # set the probability that the path continues beyond U
#     prob_continue_path = 0
#     for i in range(U, num_steps):
#         prob_continue_path += If(S[i][-1] != theta[-1], prob_path_encoding * (1-P)**(i-L), 0)
    
#     # set the total probability of the encoding
#     prob_encoding_total = prob_encoding_init * prob_path_encoding * prob_target_reached * prob_continue_path
    
#     # return the encoding
#     return And(enc_0, path_encoding, target_reached, between_L_U, prob_encoding == prob_encoding_total)
