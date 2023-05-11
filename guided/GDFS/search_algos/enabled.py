
def enabled (model, diag, nodes, curr_node, path_prob, prob_threshold):
    R = set()
    for i, r in enumerate(model.get_reactions_vector()):
        