import re
import stormpy
import stormpy.core
import __future__

def format_expression(expression):
    expression = re.sub(r'\^' , '**', expression)
    expression = re.sub(r'\&' , 'and', expression)
    expression = re.sub(r'\|' , 'or', expression)
    expression = re.sub(r'\b%s\b' % '\=', '==', expression)
    expression = re.sub(r'\b%s\b' % 'true', 'True', expression)
    return expression

def evaluate_compiled_expression(expression, var_assignments):
    return eval(expression, var_assignments)

class PRISM_Parser:
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
                    expression = format_expression(expression)
                    expression = compile(expression, '<string>', 'eval', __future__.division.compiler_flag)
                    self.reactions_vector[i][var_index] = evaluate_compiled_expression(expression,assignment)
        
        print(self.index_to_reaction_dict)
        #reaction_rate_expression is a dictionary with a list of size two for every reaction index
        #element1: compiled guard expression
        #element2: compiled rate expression
        #element3: non-compiled guard expression
        self.reaction_rate_expression = {}

        for i, c in enumerate(self.commands):
            index = self.reaction_to_index_dict[c]
            guard = ""
            flag = True
            for g in self.commands[c][0]:
                if flag:
                    guard = str(g)
                    flag = False
                else:
                    guard = guard + " and " + str(g)

            rate = ""
            flag = True
            for u in self.commands[c][1]:
                if flag:
                    rate = str(u.probability_expression)
                    flag = False
                else:
                    rate = rate + " * " + str(u.probability_expression)
            
            guard = format_expression(guard)
            guard_temp = guard
            guard = compile(guard, '<string>', 'eval', __future__.division.compiler_flag)
            rate = format_expression(rate)
            rate = compile(rate, '<string>', 'eval', __future__.division.compiler_flag)
            self.reaction_rate_expression[index] = [guard, rate, guard_temp]

        # reactions_lhs_vector contains the effect of the lhs of a 
        # reaction (not the final effect after production happens)
        self.reactions_lhs_vector = [None] * len(self.commands)
        for i, c in enumerate(self.commands):
            self.reactions_lhs_vector[i] = [None] * len(self.species_tuple)
            delimiter = 'and'
            result = [part.strip() for part in self.reaction_rate_expression[i][2].split(delimiter)]
            for r in result:
                r = re.sub(r'\(|\)', '', r)
                for j, s in enumerate(self.species_tuple):
                    temp = r[0:r.find(' ')]

                    val = 0
                    r_compiled = compile(r, '<string>', 'eval', __future__.division.compiler_flag)
                    if s == temp and '*' not in r:    
                        while (evaluate_compiled_expression(r_compiled, {s : val}) == False):
                            val = val + 1
                    if self.reactions_lhs_vector[i][j] == None:
                        self.reactions_lhs_vector[i][j] = val
                    else:
                        if val > self.reactions_lhs_vector[i][j]:
                            self.reactions_lhs_vector[i][j] = val



        

    def get_initial_state(self):
        return self.initial_state_tuple

    def get_species_tuple(self):
        return self.species_tuple

    def get_reactions_vector(self):
        return self.reactions_vector

    # def get_reactions_lhs_vector(self):
    #     return self.reactions_lhs_vector

    def get_reaction_rate(self, var_assignment_tuple, r_index):
        var_dict = {}
        for i, s in enumerate(self.species_tuple):
            var_dict[s] = var_assignment_tuple[i]
        
        guard = self.reaction_rate_expression[r_index][0]
        rate = self.reaction_rate_expression[r_index][1]

        if not evaluate_compiled_expression(guard, var_dict):
            return 0.0
        return evaluate_compiled_expression(rate, var_dict)
    
    def get_reaction_rate_constants(self):
        var_dict = {}
        for i, s in enumerate(self.species_tuple):
            var_dict[s] = 1
        
        rate_constants = []
        for r_index, r in enumerate(self.reactions_vector):
            rate = self.reaction_rate_expression[r_index][1]
            rate_constants.append(evaluate_compiled_expression(rate, var_dict))


        return rate_constants
