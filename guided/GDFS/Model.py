import stormpy
import stormpy.core
import __future__

def evaluate_expression_arithmetic(expression, variables):
    expression = expression.replace('^', '**')
    # Create a replacement dictionary from the variables dictionary
    replacement_dict = {var: str(val) for var, val in variables.items()}
    # Replace variables in the expression string using the replacement dictionary
    for var, val in replacement_dict.items():
        expression = expression.replace(var, val)
    # Evaluate the expression and return the result
    return eval(compile(expression, '<string>', 'eval', __future__.division.compiler_flag))

def evaluate_expression_logical(expression, variables):
    expression = expression.replace('&', 'and')
    expression = expression.replace('|', 'or')
    expression = expression.replace(' = ', ' == ')
    expression = expression.replace('true', 'True')
    # Create a replacement dictionary from the variables dictionary
    replacement_dict = {var: str(val) for var, val in variables.items()}
    # Replace variables in the expression string using the replacement dictionary
    for var, val in replacement_dict.items():
        expression = expression.replace(var, val)
    # Evaluate the expression and return the result
    return eval(expression)

class Model:
    def __init__(self, model_path):
        prism_program = stormpy.parse_prism_program(path=model_path, prism_compat=True)
        self.initial_state_tuple = tuple()
        vars = {}
        for m in prism_program.modules:
            for b in m.boolean_variables:
                raise Exception("|| module {} has a boolean variable - not supported ||".format(m.name))
            for v in m.integer_variables:
                vars[v.name] = v.initial_value_expression.evaluate_as_int()
        self.species_tuple = tuple(vars.keys())
        self.initial_state_tuple = tuple(vars.values())
        self.species_to_index_dict = {}
        self.index_to_species_dict = {}
        for i, s in enumerate(self.species_tuple):
            self.species_to_index_dict[s] = i
            self.index_to_species_dict[i] = s

        #commands: dictionary : {"string" : [[guards][updates]]}
        #dictionary: keys: command label | values: a list of size 2
        #first element: a list with all the guards for commands labeled the same way
        #second element: a list with all the updates for commands labeled the same way
        self.commands = {}
        for m in prism_program.modules:
            for c in m.commands:
                if len(c.updates) != 1:
                    raise Exception("|| a command has more than/less than 1 update - not supported ||")
                if len(c.action_name)==0:
                    raise Exception("|| a command doesn't have any labels - not supported ||")
                for u in c.updates:
                    if c.action_name in self.commands:
                        self.commands[c.action_name][1].append(u)
                    else:
                        self.commands[c.action_name] = [[],[u]]
                self.commands[c.action_name][0].append(c.guard_expression)
    
        self.reaction_to_index_dict = {} #index --> reaction name (label)
        self.index_to_reaction_dict = {} #reaction name --> index
        self.reactions_vector = [None] * len(self.commands)
        for i, c in enumerate(self.commands):
            self.reaction_to_index_dict[c] = i
            self.index_to_reaction_dict[i] = c
            self.reactions_vector[i] = [None] * len(self.species_tuple)
            for j, e in enumerate(self.reactions_vector[i]):
                self.reactions_vector[i][j] = 0
            for u in self.commands[c][1]:
                for a in u.assignments:
                    var_name = a.variable.name
                    var_index = self.species_to_index_dict[var_name]
                    assignment = {var_name : 0}
                    expression = str(a.expression)
                    self.reactions_vector[i][var_index] = evaluate_expression_arithmetic(expression,assignment)

    def get_initial_state(self):
        return self.initial_state_tuple

    def get_species_tuple(self):
        return self.species_tuple

    def get_reactions_vector(self):
        return self.reactions_vector

    def get_reaction_rate(self, var_assignment_tuple, r_index):
        r_label = self.index_to_reaction_dict[r_index]
        return_value = 1.0
        var_dict = {}
        for i, s in enumerate(self.species_tuple):
            var_dict[s] = var_assignment_tuple[i]
        c = self.commands[r_label]
        for g in c[0]:
            if not evaluate_expression_logical(str(g), var_dict):
                return 0.0
        for u in c[1]:
            return_value = return_value * (evaluate_expression_arithmetic(str(u.probability_expression), var_dict))
        return return_value
