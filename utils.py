from z3 import *
from Graph import graph

def exclude_path(path): 
	assignments = []
	for d in path.decls(): 
		x = Int(d.name())
		assignments.append(x == path[d])
	return Not(And(assignments))

#get smt encoding for the initial state constraint
def get_initial_state(model):
	species_vector = model.species_vector()
	initial_state = model.initial_state()
	constraints = []
	for i, s in enumerate(species_vector): 
		var_name = s + '.' + '0'
		x = Int(var_name)
		init_value = initial_state[i]
		constraints.append(x == init_value)
	return And(constraints)

def get_encoding(model, bound):
	species_vector = model.species_vector()
	constraints = []
	for j, r in enumerate(model.reactions_vector()): 
		r_constraints = []
		for i, s in enumerate(species_vector): 
			var_name_prev = s + '.' + str(bound-1)
			var_name_curr = s + '.' + str(bound)
			x = Int(var_name_prev)
			y = Int(var_name_curr)
			r_constraints.append(y == (x + model.reactions_vector()[r][2][i]))
			r_constraints.append(y >= 0)
			reaction_var_name = 'selected_reaction.' + str(bound-1) 
			selected_reaction = Int(reaction_var_name)
		r_constraints.append(selected_reaction == (j+1))
		constraints.append(And(r_constraints))
	return simplify(Or(constraints))


def construct_path_o(graph, model, model_name, prism, csl_prop, cp_bound, count_, prob_, prob_bound, property_var, property_val, mc_step):
	count = count_
	prob = prob_
	curr_bound = 1
	c_count = 0
	#initial state of a path can be any node in the graph
	init_const = []
	for n in graph.node_list:
		var_values = []
		if n.var_dict[property_var] == property_val: 
			continue
		for s in n.var_dict: 
			x = Int (s + '.0')
			temp_const = (x == n.var_dict[s])
			var_values.append(temp_const)
		init_const.append(And(var_values))
	init_const = Or(init_const)

	while curr_bound<cp_bound and c_count<50:
		flag = False
		solver = Solver()
		solver.add(init_const)
		solver.add(exclude_graph(graph, curr_bound))
		for j in range(1, curr_bound):
			solver.add(loop_constraint(model, j))
			solver.add(get_encoding(model, j))
		solver.add(get_encoding(model, curr_bound))

		#target states can be any state on the graph or a new
		#target state
		target_const_temp = []
		for n in graph.node_list:
			var_values = []
			for s in n.var_dict:
				x = Int (s + '.' + str(curr_bound))
				temp_const = (x == n.var_dict[s])
				var_values.append(temp_const)
			target_const_temp.append(And(var_values))
		target_const_temp = Or(target_const_temp)
		x = Int(property_var + '.' + str(curr_bound))
		property_constraint = (x==property_val)
		target_const = Or(target_const_temp, property_constraint)

		while(solver.check(target_const)==sat):
			c_count = c_count + 1
			if c_count > 50: 
				break
			flag = True
			path = solver.model()
			graph.add_path_o(path)
			solver.add(exclude_path(path))
			count = count + 1
			# if (count%mc_step)==0: 
			# 	prob = graph.model_check(model, model_name, prism, csl_prop)
			# 	#print(prob)
			# 	# print('total number of paths so far: ' + str(count))
			# 	# print('probability= ' + str(prob))
			# 	# print('+++++')
			# if prob>=prob_bound:
			# 	return count, prob, True

		curr_bound = curr_bound + 1
		if flag: 
			curr_bound = 1

	return count, prob, False


def get_encoding_o(model, bound, prob_growth, reaction_sets, current):
	species_vector = model.species_vector()
	constraints = []
	for j, r in enumerate(model.reactions_vector()): 
		r_constraints = []
		for i, s in enumerate(species_vector): 
			var_name_prev = s + '.' + str(bound-1)
			var_name_curr = s + '.' + str(bound)
			x = Int(var_name_prev)
			y = Int(var_name_curr)
			r_constraints.append(y == (x + model.reactions_vector()[r][2][i]))
			r_constraints.append(y >= 0)
			reaction_var_name = 'selected_reaction.' + str(bound-1) 
			selected_reaction = Int(reaction_var_name)
		r_constraints.append(selected_reaction == (j+1))
		constraints.append(And(r_constraints))
	return simplify(Or(constraints))