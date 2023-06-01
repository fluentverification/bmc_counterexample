import itertools
from z3 import Solver, Int, And, sat
from Parser import Parser
import time
import json
import subprocess
import copy

#
def CEX_GEN(json_data):
    start_time = time.time()
    
    #parsing parameters from json elements
    jani_path = json_data['jani_path']
    model_name = json_data['model_name']
    model_path = json_data['model_path']
    K = int(json_data['starting_bound'])
    csl_prop_lb = json_data['csl_property']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    max_comb = int(json_data['max_combination'])
    storm_bin = json_data["storm_binary"]
    #

    #parse the model into Parser() object
    model = Parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    #
    
    #generate subsets (keys for the min_max dictionary)
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    subsets = get_subsets(index_tuple, 1, max_comb)
    #

    print("\nStarting...")
    while (True):
        print("Bound = " + str(K))
        before = time.time()
        min_max_dict = get_min_max(model, K, target_index, target_value, subsets)
        print("Generating min_max dictonary took " + str(time.time() - before) + "seconds.")
        JSON_Parser(model, model_name, K, jani_path, min_max_dict)
        print("Calling Storm to calculate the probability... \n\n")
        
        #running storm on the produced output
        stdout_result = subprocess.run([storm_bin, "--jani", "./results/" + model_name + "/bounds/" + model_name + "_" + str(K) + ".jani", '--prop', csl_prop_lb], stdout=subprocess.PIPE)
        stdout_result = stdout_result.stdout.decode('utf-8')
        print(stdout_result)
        #
        
        K = K + 1
        print("\n \nRunning time: " + str(time.time() - start_time) + " seconds")
        print("=" * 50)
        # break
#

#
def JSON_Parser(model, model_name, K, jani_path, min_max_dict):

    try:
        with open(jani_path, "r") as json_file:
            parsed_json = json.load(json_file)
    except FileNotFoundError:
        print(f"File '{jani_path}' not found.")
    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON: {e}")


    species_tuple = model.get_species_tuple()
    bounds_tuple = ([None, None], )
    for i in range(1,len(model.get_species_tuple())):
        bounds_tuple = bounds_tuple + ([None, None],)
    
    for e in min_max_dict:
        if len(e) == 1:
            for i, s in enumerate(species_tuple):
                if i==e[0]:
                    bounds_tuple[i][0] = min_max_dict[e][0]
                    bounds_tuple[i][1] = min_max_dict[e][1]
    

    global_variables = parsed_json["variables"]
    automata = parsed_json["automata"]

    #adding lower-bound and upper-bound to global variables of the model
    for gv in global_variables:
        for i, s in enumerate(species_tuple):
            if gv["name"] == s:
                type_ = {"base" : "int", 
                "kind" : "bounded", 
                "lower-bound" : bounds_tuple[i][0],
                "upper-bound" : bounds_tuple[i][1]
                }
                gv["type"] = type_
    #

    #adding lower-bound and upper-bound to local variables of an automaton (module)
    for automaton in automata:
        if "variables" in automaton:
            local_variables = automaton["variables"]
            for lv in local_variables:
                for i, s in enumerate(species_tuple):
                    if lv["name"] == s:
                        type_ = {"base" : "int", 
                        "kind" : "bounded", 
                        "lower-bound" : bounds_tuple[i][0],
                        "upper-bound" : bounds_tuple[i][1]
                        }
                        lv["type"] = type_
    #

    #adding a new sink variable to the model
    global_variables = parsed_json["variables"]
    sink_variable = {
            "initial-value": 0,
            "name": "sink_var",
            "type": {
                "base": "int",
                "kind": "bounded",
                "lower-bound": 0,
                "upper-bound": 1
            }
        }
    global_variables.append(sink_variable)
    #
    
    #generating the guards that guarantee reducing variables does not result in variables value 
    #dropping below its lower-bound
    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            destinations = edge["destinations"]
            for destination in destinations:
                if "assignments" in destination:
                    assignments = destination["assignments"]
                    lower_bound_guards = []
                    for assignment in assignments:
                        for i, s in enumerate(species_tuple):
                            if assignment["ref"] == s:
                                value = assignment["value"]
                                if "op" not in value:
                                    raise Exception("unsupported value assignment")
                                if value["op"] == "-":
                                    rhs = value["right"]
                                    gt_guard = {"left" : s, "op" : ">", "right" : bounds_tuple[i][0] + rhs - 1}
                                    lower_bound_guards.append(gt_guard)
            if "guard" not in edge:
                raise Exception("an edge does not have a guard")
            guard = edge["guard"]
            flag = False
            for g in lower_bound_guards:
                if not flag:
                    guard["exp"] = {"left" : guard["exp"], "op": "∧", "right" : g}
                else:
                    guard["exp"] = {"left" : guard["exp"]["left"], "op" : "∧", "right" : {"left": guard["exp"]["right"], "op" : "∧", "right" : g}}
                flag = True
            if flag:
                edge["comment"] = "modified"
                guard["comment"] = "modified! old comment: " + guard["comment"]
    #

    #generating the guards that guarantee increasing variables does not result in variables value 
    #going above its upper-bound
    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            destinations = edge["destinations"]
            for destination in destinations:
                if "assignments" in destination:
                    assignments = destination["assignments"]
                    upper_bound_guards = []
                    for assignment in assignments:
                        for i, s in enumerate(species_tuple):
                            if assignment["ref"] == s:
                                value = assignment["value"]
                                if "op" not in value:
                                    raise Exception("unsupported value assignment")
                                if value["op"] == "+":
                                    rhs = value["right"]
                                    lt_guard = {"left" : s, "op" : "<", "right" : bounds_tuple[i][1] - rhs + 1}
                                    upper_bound_guards.append(lt_guard)
            if "guard" not in edge:
                raise Exception("an edge does not have a guard")
            guard = edge["guard"]
            flag = False
            for g in upper_bound_guards:
                if (not flag) and ("modified" not in guard["comment"]):
                    guard["exp"] = {"left" : guard["exp"], "op": "∧", "right" : g}
                else: 
                    guard["exp"] = {"left" : guard["exp"]["left"], "op" : "∧", "right" : {"left": guard["exp"]["right"], "op" : "∧", "right" : g}}
                flag = True
            if flag:
                edge["comment"] = "modified"
                guard["comment"] = "modified! old comment: " + guard["comment"]
    #

    #expanding the guards for all edges with min_max constraints
    min_max_exp = {}
    for e in min_max_dict:
        lhs = {}
        if len(e) == 1:
            lhs = {"left" : species_tuple[e[0]], "op" : "+", "right" : 0} 
        else:
            flag = True
            for i in e:
                if flag:
                    lhs = {"left" : species_tuple[i], "op" : "+", "right" : 0}
                    flag = False
                else:
                    lhs = {"left" : lhs, "op" : "+", "right" : species_tuple[i]}
        if len(min_max_exp) == 0:
            min_max_exp = {"left" : {"left" : lhs, "op": ">", "right" : min_max_dict[e][0] - 1}, 
            "op": "∧", 
            "right": {"left" : lhs, "op": "<", "right" : min_max_dict[e][1] + 1}
            }
        else:
            min_max_exp = {"left" : min_max_exp, "op" : "∧", 
            "right" : { "left" : {"left" : lhs, "op": ">", "right" : min_max_dict[e][0] - 1}, 
            "op": "∧", 
            "right": {"left" : lhs, "op": "<", "right" : min_max_dict[e][1] + 1}
            }
            }

    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            if "guard" not in edge:
                raise Exception("an edge does not have a guard")
            guard = edge["guard"]
            if ("modified" not in guard["comment"]):
                guard["exp"] = {"left" : guard["exp"], "op": "∧", "right" : min_max_exp}
            else: 
                guard["exp"] = {"left" : guard["exp"]["left"], "op" : "∧", "right" : {"left": guard["exp"]["right"], "op" : "∧", "right" : min_max_exp}}
            guard["comment"] = "modified! old comment: " + guard["comment"]
    #
                        
    #Iterating over all the edges in automata and generating a new edge 
    #going to sink state for each modified edge in order to maintain
    #semantics of the original model
    
    ##1-generating the sink state:
    sink_assignments = []
    sink_assignments.append({"ref" : "sink_var", "value" : 1})
    for i, s in enumerate(species_tuple):
        sink_assignments.append({"ref" : s, "value" : bounds_tuple[i][0]})
    ##

    for automaton in automata:
        edges = automaton["edges"]
        new_edges = []
        for edge in edges:
            if "comment" not in edge:
                continue
            if edge["comment"] != "modified":
                continue
            new_edge = copy.deepcopy(edge)
            new_edge["comment"] = "semantics edge"
            guard = new_edge["guard"]
            guard["comment"] = "semantics guard"
            guard["exp"] = {"left" : guard["exp"]["left"], "op": "∧", "right" : {"op" : "¬", "exp" : guard["exp"]["right"]}}
            destinations = new_edge["destinations"]
            for destination in destinations:
                if "assignments" in destination:
                    destination["assignments"] = sink_assignments
            new_edges.append(new_edge)
        for new_edge in new_edges:
            edges.append(new_edge)
    #
    
    #Exporting the modified model for a bound into a jani file
    file_path = "./results/" + model_name + "/bounds/" + model_name + "_" + str(K) + ".jani"
    try:
        with open(file_path, "w") as json_file:
            json.dump(parsed_json, json_file, indent=4, ensure_ascii=False)
        print(f"JSON data saved to '{file_path}' successfully.")
    except IOError:
        print(f"Unable to save JSON data to '{file_path}'.")
    #
#





############ Functions Generating Min_Max via SMT-Solving ############
#returns all the subsets of a tuple with cardinality in range [L, U]
#returns a list of tuples where each tuple is a subset
def get_subsets (tuple, L, U):
    if U>len(tuple):
        U = len(tuple)
    return list(itertools.chain.from_iterable(itertools.combinations(tuple, r) for r in range(L, U+1)))
#

#generating a dictionary. 
#Keys: combinations of variables:(s1,), (s1,s2), ...
#Values: a list with l[0] the minimum value and l[1] the maximum value for that key
def get_min_max(model, bound, target_index, target_value, subsets):
    min_max = {}
    index_tuple = (0,)
    for i in range(1,len(model.get_species_tuple())):
        index_tuple = index_tuple + (i,)
    for s in subsets:
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
        if r[target_index]!=0:
            vars.append([Int("n_i_" + str(i)), r[target_index]])
            vars.append([Int("n_t_" + str(i)), r[target_index]])
    sum = model.get_initial_state()[target_index]
    for i in vars:
        sum = sum + i[0]*i[1]
    constraints.append(sum==target_value)
    
    #fourth constraint (species population>=0)
    for i, s in enumerate(model.get_species_tuple()):
        vars_i = []
        vars = []
        for j, r in enumerate(model.get_reactions_vector()):
            if r[i]!=0:
                vars.append([Int("n_i_"+str(j)), r[i]])
                vars.append([Int("n_t_"+str(j)), r[i]])
                vars_i.append([Int("n_i_"+str(j)), r[i]])
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
                    max_value = max_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][1] = max_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
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
                    min_value = min_value + (model.get_reactions_vector()[v[0]][j] * v[1])
            min_max[e][0] = min_value
            sum = 0
            for j in e:
                sum = sum + model.get_initial_state()[j]
            for v in vars:
                for j in e:
                    sum = sum + (Int(v[2]) * model.get_reactions_vector()[v[0]][j])
            solver.add(sum<min_value)
        solver.pop()
    return min_max
#