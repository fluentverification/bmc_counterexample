#! /home/ubu/projects/probmc/storm-1.9/env/bin/python3
import stormpy

class prism_program:
    def __init__(self, prism_program_path):
        self.pp = stormpy.parse_prism_program(prism_program_path, prism_compat = True)

        #collecting all the variables in the program{
        self.variables = []
        for v in self.pp.global_integer_variables:
            self.variables.append(v.expression_variable)
        for v in self.pp.global_boolean_variables:
            self.variables.append(v.expression_variable)
        for m in self.pp.modules:
            for v in m.integer_variables:
                self.variables.append(v.expression_variable)
            for v in m.boolean_variables:
                self.variables.append(v.epxression_variable)
        #}


path = "./stormpy_test/enzym_unb.sm"        
pp = prism_program(path) 
for v in pp.variables:
    print (v.name)
    print(type(v))
    if v.has_integer_type():
        print("integer_variable")
    elif v.has_boolean_type():
        print("boolean_variable")
    else:
        raise Exception("unsupported variable type in the prism program (other than int/bool)")
    print(v.initial_value_expression)
    print("---")
quit()
prism_program = stormpy.parse_prism_program("./enzym_unb.sm", prism_compat = True)
#print("#"*30)
#print(prism_program.nr_modules)
#print(prism_program.modules)
#print(prism_program.has_initial_states_expression)
#print(prism_program.variables)
for v in prism_program.variables:
    #print(v.name)
    pass
for c in prism_program.constants:
    #print(c.name)
    pass
print("="*10)

#collecting all the variables in the program{
variables = []
for v in prism_program.global_integer_variables:
    variables.append(v)
for v in prism_program.global_boolean_variables:
    variables.append(v)
for m in prism_program.modules:
    for v in m.integer_variables:
        variables.append(v)
    for v in m.boolean_variables:
        variables.append(v)
for v in variables:
    print(str(type(v)) + " : "  + str(v.initial_value_expression))
#}
print("then:")
#getting the initial state of the program as a vector of bools and ints{
initial_state_stormpy = [None] * len(variables)
for i, v in enumerate(variables):
    initial_state_stormpy[i] = v.initial_value_expression
if None in initial_state_stormpy:
    raise Exception("initial state was not retrieved")
initial_state = [None] * len(variables)
for i, v in enumerate(initial_state_stormpy):
    if v.has_integer_type():
        initial_state[i] = v.evaluate_as_int()
    elif v.has_boolean_type():
        initial_state[i] = v.evaluate_as_bool()
    else:
        raise Exception("one of the variables in the model has a type other than bool/int")
print(initial_state)
#}
print("then:")
#getting the commands in the program{
commands = {}
count = 1
for m in prism_program.modules:
    for c in m.commands:
        if c.is_labeled:
            if c.action_name in commands.keys():
                commands[c.action_name].append(c)
            else:
                commands[c.action_name] = []
                commands[c.action_name].append(c)
        else:
            commands["temp_label_" + str(count)] = c
            count = count + 1
print(commands)
#}
#next_states: curr_state, commands -> [is_enabled, next_state]^commands
def next_states(curr_state, commands):
    to_return = [None] * len(commands)
    for i, command in enumerate(commands):
        if is_enabled (command, curr_state):
            to_return[i] = [True, X]
        else:
            to_return[i] = [False, None]
#}
#is_enabled: (curr_state, command) -> {True, False} {
def is_enabled(curr_state, command, variables):
    substitute_dict = {}
    for i, v in enumerate(variables):
        substitute_dict[v.expression_variable] = stormpy.ExpressionManager().create_integer(curr_state[i])
    to_return = True
    print(substitute_dict)
    for c in commands[command]:
        to_return = to_return and c.guard_expression.substitute(substitute_dict).evaluate_as_bool()
    return to_return
#}

#module = prism_program.get_module(name=)
print("#"*30)
for command in commands:
    print (command)
    print (is_enabled(initial_state, command, variables))
    print("***")
quit()

# formula_str = "P=? [F<=100 Gbg>=50]"
# properties = stormpy.parse_properties(formula_str, prism_program)
# model = stormpy.build_model(prism_program, properties)
# print("Number of states: {}".format(model.nr_states))
# print("Number of transitions: {}".format(model.nr_transitions))
# 
# result = stormpy.model_checking(model, properties[0])
# initial_state = model.initial_states[0]
# print(result.at(initial_state))
# 
#
