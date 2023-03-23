def get_total_outgoing_rate(var_values, model):
    rate = 0
    for i, r in enumerate(model.reactions_vector()):
        comb = 1
        for j, coefficient in enumerate(r[0]):
            if coefficient>0:
                for c in range(coefficient):
                    comb = comb * var_values[j]
        rate = rate + model.reaction_rates()[i] * comb
    return rate

def get_reaction_rate(var_values, model, reaction):
    rate = model.reaction_rates()[reaction]
    for i, coefficient in enumerate(model.reactions_vector()[reaction][0]):
        if coefficient>0:
            for c in range(coefficient):
                rate = rate * var_values[i]
    return rate

class set_(set):
    def __contains__(self, value):
        # Override the __contains__ method
        return any(value.equals(x) for x in self)