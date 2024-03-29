import copy
import json
import math


def JANI_Parser(model, model_name, jani_model, min_max_dict):
    try:
        with open(jani_model, "r") as json_file:
            parsed_json = json.load(json_file)
    except FileNotFoundError:
        print(f"File '{jani_model}' not found.")
    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON: {e}")
    
    #1-add the bounds for each single species to the variables
    parsed_json = add_bounds(parsed_json=parsed_json, model=model, min_max_dict=min_max_dict)

    #2-population of species should remain within provided guards
    parsed_json = population_guard(parsed_json=parsed_json, model=model, min_max_dict=min_max_dict)

    #3-adding the constraints on population of multiple species
    parsed_json = multiple_species_guard(parsed_json=parsed_json, model=model, min_max_dict=min_max_dict)

    #4-add a sink variable and define a sink state representing the probability of the rest of the state space
    parsed_json, sink_state = sink_assignment(parsed_json=parsed_json)
    
    #5-reactions can only fire from non-sink states. Adding sink_var==0 to the guard of all reactions.
    parsed_json = sink_guard(parsed_json=parsed_json)
    
    #6-For each reaction with augmented guard, add a new reaction where:
    #new guard = (original guard) ^ (! additional guard)
    #destination = sink state
    #rate = old rate
    parsed_json = semantic_guard(parsed_json=parsed_json, sink_assignment=sink_state)
    
    #Exporting the modified model into a jani fil with sink states
    file_path = "./tmp/" + model_name + "_sink" + ".jani"
    try:
        with open(file_path, "w") as json_file:
            json.dump(parsed_json, json_file, indent=4, ensure_ascii=False)
        print(f"JSON data saved to '{file_path}' successfully.")
    except IOError:
        print(f"Unable to save JSON data to '{file_path}'.")
    #

#
        
def add_bounds(parsed_json, model, min_max_dict):
    #getting the minimum and maximum values for species in all bounds
    species_tuple = model.get_species_tuple()
    var_bounds_tuple = ([None, None], )
    for i in range(1,len(model.get_species_tuple())):
        var_bounds_tuple = var_bounds_tuple + ([None, None],) 
    for k, v in min_max_dict.items():
        for key, value in v.items():
            if len(key) == 1:
                for i, s in enumerate(species_tuple):
                    if i==key[0]:
                        if var_bounds_tuple[i][0] == None:
                            var_bounds_tuple[i][0] = value[0]
                        elif value[0] < var_bounds_tuple[i][0]:
                            var_bounds_tuple[i][0] = value[0]
                        if var_bounds_tuple[i][1] == None:
                            var_bounds_tuple[i][1] = value[1]
                        elif value[1] > var_bounds_tuple[i][1]:
                            var_bounds_tuple[i][1] = value[1]
    #
    
    global_variables = parsed_json["variables"]
    automata = parsed_json["automata"]
    
    #adding lower-bound and upper-bound to global variables of the model
    for gv in global_variables:
        for i, s in enumerate(species_tuple):
            if gv["name"] == s:
                type_ = {"base" : "int", 
                "kind" : "bounded", 
                "lower-bound" : var_bounds_tuple[i][0],
                "upper-bound" : var_bounds_tuple[i][1]
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
                        "lower-bound" : var_bounds_tuple[i][0],
                        "upper-bound" : var_bounds_tuple[i][1]
                        }
                        lv["type"] = type_
    #
    
    #after all the bounds were added, if there still exists an unbounded variable,
    #choose its initial value as both the lower_bound and the upper_bound
    #(these are constants defined as variables)
    for gv in global_variables:
        if gv["type"] == "int":
            type_ = {"base" : "int", 
                    "kind" : "bounded", 
                    "lower-bound" : gv["initial-value"],
                    "upper-bound" : gv["initial-value"]
                    }
            gv["type"] = type_
    
    for automaton in automata:
        if "variables" in automaton:
            local_variables = automaton["variables"]
            for lv in local_variables:
                if lv["type"] == "int":
                    type_ = {"base" : "int", 
                            "kind" : "bounded", 
                            "lower-bound" : lv["initial-value"],
                            "upper-bound" : lv["initial-value"]
                            }
                    lv["type"] = type_
    #
    
    return parsed_json
#

def population_guard(parsed_json, model, min_max_dict):
    automata = parsed_json["automata"]
    species_tuple = model.get_species_tuple()
    #getting the minimum and maximum values for species in all bounds
    bounds_tuple = ([None, None], )
    for i in range(1,len(model.get_species_tuple())):
        bounds_tuple = bounds_tuple + ([None, None],) 
    for k, v in min_max_dict.items():
        for key, value in v.items():
            if len(key) == 1:
                for i, s in enumerate(species_tuple):
                    if i==key[0]:
                        if bounds_tuple[i][0] == None:
                            bounds_tuple[i][0] = value[0]
                        elif value[0] < bounds_tuple[i][0]:
                            bounds_tuple[i][0] = value[0]
                        if bounds_tuple[i][1] == None:
                            bounds_tuple[i][1] = value[1]
                        elif value[1] > bounds_tuple[i][1]:
                            bounds_tuple[i][1] = value[1]
    #

    #generating the guards that guarantee reducing variables does not result in variables value 
    #dropping below its lower-bound
    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            if "guard" not in edge:
                # raise Exception("an edge does not have a guard")
                guard = {'comment': 'generated empty guard', 'exp': True}
                edge["guard"] = guard
            else:
                guard = edge["guard"]

            destinations = edge["destinations"]
            for destination in destinations:
                if "assignments" in destination:
                    assignments = destination["assignments"]
                    for assignment in assignments:
                        for i, s in enumerate(species_tuple):
                            if assignment["ref"] == s:
                                value = assignment["value"]
                                if "op" not in value:
                                    raise Exception("unsupported value assignment")
                                if value["op"] == "-":
                                    edge["comment"] = "modified"
                                    rhs = value["right"]
                                    gt_guard = {"left" : bounds_tuple[i][0] + rhs, "op" : "≤", "right" : s}
                                    if ("modified" not in guard["comment"]):
                                        guard["exp"] = {"left" : guard["exp"], "op": "∧", "right" : gt_guard}
                                        guard["comment"] = "modified"
                                    else:
                                        guard["exp"] = {"left" : guard["exp"]["left"], 
                                                        "op": "∧", 
                                                        "right" : {
                                                            "left" : guard["exp"]["right"],
                                                            "op" : "∧",
                                                            "right" : gt_guard
                                                        }
                                                        }
                                        
                                if value["op"] == "+":
                                    edge["comment"] = "modified"
                                    rhs = value["right"]
                                    lt_guard = {"left" : s , "op" : "≤", "right" : bounds_tuple[i][1] - rhs}
                                    if ("modified" not in guard["comment"]):
                                        guard["exp"] = {"left" : guard["exp"], "op": "∧", "right" : lt_guard}
                                        guard["comment"] = "modified"
                                    else:
                                        guard["exp"] = {"left" : guard["exp"]["left"], 
                                                        "op": "∧", 
                                                        "right" : {
                                                            "left" : guard["exp"]["right"],
                                                            "op" : "∧",
                                                            "right" : lt_guard
                                                        }
                                                        }
                                        

    #
    return parsed_json
#

def multiple_species_guard(parsed_json, model, min_max_dict):
    species_tuple = model.get_species_tuple()
    min_max_exp = {}
    for k, e in min_max_dict.items():
        min_max_bound = {}
        for key, element in e.items():
            if len(key) == 1:
                continue
            flag = True
            for i in key:
                if flag:
                    lhs = species_tuple[i]
                    flag = False
                else:
                    lhs = {"left" : lhs, "op" : "+", "right" : species_tuple[i]}
            if len(min_max_bound) == 0 :
                min_max_bound = {"left" : {"left" : element[0], 
                                           "op": "≤", 
                                           "right" : lhs
                                           }, 
                                "op": "∧", 
                                "right": {"left" : lhs, 
                                          "op": "≤", 
                                          "right" : element[1]}
                                }
            else:
                min_max_bound = {"left" : min_max_bound,
                                 "op" : "∧",
                                 "right" : { "left" : { "left" : element[0], 
                                                        "op": "≤", 
                                                        "right" : lhs
                                                        }, 
                                            "op": "∧", 
                                            "right": {"left" : lhs, 
                                                        "op": "≤", 
                                                        "right" : element[1]
                                                    }    
                                            }
                                 }

        if len(min_max_exp) == 0:
            min_max_exp = min_max_bound
        else:
            min_max_exp = {"left" : min_max_exp,
                           "op" : "∨",
                           "right" : min_max_bound
                           }

    
        
    automata = parsed_json["automata"]
    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            guard = edge["guard"]
            if ("modified" in guard["comment"]):
                guard["exp"] = {"left" : guard["exp"]["left"],
                                "op": "∧", 
                                "right" : {"left" : guard["exp"]["right"],
                                           "op": "∧",
                                           "right" : min_max_exp}
                                           }
    #
                                        

    #
    return parsed_json
#

def sink_assignment(parsed_json):
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
    #generating the sink state:
    sink_assignments = []
    sink_assignments.append({"ref" : "sink_var", "value" : 1})
    for gv in global_variables:
        if (gv["name"] != "sink_var"):
            value = gv["initial-value"]
            sink_assignments.append({"ref" : gv["name"], "value" : value})

    automata = parsed_json["automata"]
    for automaton in automata:
        if "variables" in automaton:
            local_variables = automaton["variables"]
            for lv in local_variables:
                if (lv["name"]!= "sink_var"):
                    value = lv["initial-value"]
                    sink_assignments.append({"ref" : lv["name"], "value" : value})
    #
    
    return parsed_json, sink_assignments
#

def sink_guard(parsed_json):
    automata = parsed_json["automata"]
    for automaton in automata:
        edges = automaton["edges"]
        for edge in edges:
            guard = edge["guard"]
            if "modified" not in guard["comment"]:
                guard["exp"] = {"left" : guard["exp"], 
                                "op" : "∧", 
                                "right" : {"left" : "sink_var",
                                                        "op" : "=",
                                                        "right" : 0
                                                        }
                                }
            elif guard["exp"] == True:
                guard["exp"] = {"left" : "sink_var", "op" : "=", "right" : 0}
            else:
                guard["exp"] = {"left" : { "left" : guard["exp"]["left"],
                                            "op" : "∧",
                                            "right" : {"left" : "sink_var",
                                                        "op" : "=",
                                                        "right" : 0
                                                        }
                                            },
                                "op" : "∧",
                                "right": guard["exp"]["right"]}
    return parsed_json
#

def semantic_guard(parsed_json, sink_assignment):
    automata = parsed_json["automata"]
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
            del new_edge["action"]
            guard = new_edge["guard"]
            guard["comment"] = "semantics guard"
            guard["exp"] = {"left" : guard["exp"]["left"], "op": "∧", "right" : {"op" : "¬", "exp" : guard["exp"]["right"]}}
            destinations = new_edge["destinations"]
            for destination in destinations:
                destination["assignments"] = sink_assignment
            new_edges.append(new_edge)
        for new_edge in new_edges:
            edges.append(new_edge)
    return parsed_json
#

