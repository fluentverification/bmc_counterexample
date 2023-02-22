from z3 import *

def check_satisfiability():
    # create z3 integer variables for x, y, and z
    x, y, z = Ints('x y z')
    
    # create a solver instance
    solver = Solver()
    
    # add the constraint that x^2 + y^2 = z^2
    solver.add(x**2 + y**2 == z**2)
    
    # check for satisfiability
    if solver.check() == sat:
        # if there is a satisfiable assignment, print the values of x, y, and z
        model = solver.model()
        print(f"x = {model[x]}, y = {model[y]}, z = {model[z]}")
    else:
        print("No satisfiable assignment exists.")

from z3 import *

def check_satisfiability_nonzero():
    # create z3 integer variables for x, y, and z
    x, y, z = Ints('x y z')
    
    # create a solver instance
    solver = Solver()
    
    # add the constraint that x^2 + y^2 = z^2
    solver.add(x**2 + y**2 == z**2)
    
    # add constraints to ensure that x, y, and z are non-zero
    solver.add(x != 0)
    solver.add(y != 0)
    solver.add(z != 0)
    
    # check for satisfiability
    if solver.check() == sat:
        # if there is a satisfiable assignment, print the values of x, y, and z
        model = solver.model()
        print(f"x = {model[x]}, y = {model[y]}, z = {model[z]}")
    else:
        print("No satisfiable assignment exists.")


for i in range(1, 5):
    print(i)

S = [ [Bool(f"S_{i}^{j}") for j in range(3)] for i in range(4) ]
print (S[0][0])

def get_encoding(i):
    # implement the encoding for taking allowable transitions until it reaches bound i
    pass

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
