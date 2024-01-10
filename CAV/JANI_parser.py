import copy
import json
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
                # raise Exception("an edge does not have a guard")
                guard = {'comment': 'generated empty guard', 'exp': True}
                edge["guard"] = guard
            else:
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
                # raise Exception("an edge does not have a guard")
                guard = {'comment': 'generated empty guard', 'exp': True}
                edge["guard"] = guard
            else:
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
                # raise Exception("an edge does not have a guard")
                guard = {'comment': 'generated empty guard', 'exp': True}
                edge["guard"] = guard
            else:
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
        sink_assignments.append({"ref" : s, "value" : bounds_tuple[i][0]-1})
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
