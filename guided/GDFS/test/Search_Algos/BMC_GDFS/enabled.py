from Utils.misc import get_total_outgoing_rate, get_reaction_rate
from math import log

# Input: diag, current node, current path probability, probability threshold
# returns a list of elements [e1, e2, e3, ...]
# each e_i is a tuple (index, flag, prob) where 'index' is the index of the 
# enabled reaction and 'flag' indicates whether the transition with that 
# reaction is already in the diag and 'prob' is the path probability ending
# in the destination of the transition. The list is sorted in ascending order
# with respect to the transition probability in the embedded-DTMC, enabling
# calling pop() on the returned list to get the reaction with highest transition
# probability
def enabled(model, diag, curr_node, path_prob, prob_thresh):
    R = []
    total_rate = get_total_outgoing_rate(curr_node.var_values, model)

    for i, r in enumerate(model.get_reactions_vector()):
        rate = get_reaction_rate(curr_node.var_values, model, i) 
        
        if rate > 0:
            p = path_prob + log(rate/total_rate)
            if p > prob_thresh:
                dst_var_values = [None] * len(model.get_species_tuple())
                for j, coefficient in enumerate(r):
                    dst_var_values[j] = curr_node.var_values[j] + coefficient
                dst_var_values = tuple(dst_var_values)

                edge_tuple = curr_node.var_values + dst_var_values + (i,)

                if diag.get_edge(edge_tuple)[0]:
                    R.append((i, False, p))
                else:
                    R.append((i, True, p))

    R =  sorted(R, key=lambda x: x[2])
    return R