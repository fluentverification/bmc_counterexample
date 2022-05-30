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

#get smt encoding for a bound
def get_encoding_2(model, bound):
	species_vector = model.species_vector()
	constraints = []
	for j, r in enumerate(model.reactions_vector()): 
		if j == 4 or j == 7 or j==2:
			r_constraints = []
			for i, s in enumerate(species_vector): 
				var_name_prev = s + '.' + str(bound-1)
				var_name_curr = s + '.' + str(bound)
				x = Int(var_name_prev)
				y = Int(var_name_curr)
				
				# if bound < 40:
				# 	r_constraints.append(y == (x + model.reactions_vector()[8][2][i]))
				# 	r_constraints.append(y >= 0)
				# 	reaction_var_name = 'selected_reaction.' + str(bound-1) 
				# 	selected_reaction = Int(reaction_var_name)
				# 	r_constraints.append(selected_reaction == 8)
				# else:
				r_constraints.append(y == (x + model.reactions_vector()[r][2][i]))

				#r_constraints.append(And(y >= 0, y<=60))
				r_constraints.append(y >= 0)
				reaction_var_name = 'selected_reaction.' + str(bound-1) 
				selected_reaction = Int(reaction_var_name)
				r_constraints.append(selected_reaction == (j+1))

			constraints.append(And(r_constraints))
	return simplify(Or(constraints))

#get smt encoding for a bound
def get_encoding_3(model, bound):
	species_vector = model.species_vector()
	constraints = []
	for j, r in enumerate(model.reactions_vector()): 
		if j == 4 or j == 7:
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

#get smt encoding for the initial state constraint
def get_inital_state_rate(model):
	species_vector = model.species_vector()
	initial_state = model.initial_state()
	constraints = []
	for i, s in enumerate(species_vector): 
		var_name = s + '.' + '0'
		x = Int(var_name)
		init_value = initial_state[i]
		constraints.append(x == init_value)
	rate = Real('rate.0')
	constraints.append(rate == 1.0)
	return And(constraints)


def get_encoding_rate(model, bound):
	species_vector = model.species_vector()
	constraints = []
	
	for j, r in enumerate(model.reactions_vector()): 
		r_constraints = []
		reactants = []
		for i, s in enumerate(species_vector): 
			if model.reactions_vector()[r][0][i] > 0: 
				x = Int(s + '.' + str(bound))
				reactants.append(x)
			var_name_prev = s + '.' + str(bound-1)
			var_name_curr = s + '.' + str(bound)
			x = Int(var_name_prev)
			y = Int(var_name_curr)
			r_constraints.append(y == (x + model.reactions_vector()[r][2][i]))
			r_constraints.append(y >= 0)
			reaction_var_name = 'selected_reaction.' + str(bound-1) 
			selected_reaction = Int(reaction_var_name)
		r_constraints.append(selected_reaction == (j+1))
		
		temp_rate = model.reaction_rates()[j] / 124.0
		temp_rate = RealVal(temp_rate)

		for s in reactants:
			temp_rate = temp_rate * s
		rate_curr = Real('rate.' + str(bound))
		rate_prev = Real('rate.' + str(bound-1))
		rate_constraint = (rate_curr == (rate_prev*temp_rate))
		r_constraints.append(rate_constraint)
		constraints.append(And(r_constraints))
	return simplify(Or(constraints))

def get_encoding_rate2(model, bound):
	species_vector = model.species_vector()
	constraints = []
	
	for j, r in enumerate(model.reactions_vector()): 
		if j == 4 or j == 7:

			r_constraints = []
			reactants = []
			for i, s in enumerate(species_vector): 
				if model.reactions_vector()[r][0][i] > 0: 
					x = Int(s + '.' + str(bound))
					reactants.append(x)
				var_name_prev = s + '.' + str(bound-1)
				var_name_curr = s + '.' + str(bound)
				x = Int(var_name_prev)
				y = Int(var_name_curr)
				r_constraints.append(y == (x + model.reactions_vector()[r][2][i]))
				r_constraints.append(y >= 0)
				reaction_var_name = 'selected_reaction.' + str(bound-1) 
				selected_reaction = Int(reaction_var_name)
			r_constraints.append(selected_reaction == (j+1))
			
			temp_rate = model.reaction_rates()[j] / 124.0
			temp_rate = RealVal(temp_rate)

			for s in reactants:
				temp_rate = temp_rate * s
			rate_curr = Real('rate.' + str(bound))
			rate_prev = Real('rate.' + str(bound-1))
			rate_constraint = (rate_curr == (rate_prev*temp_rate))
			r_constraints.append(rate_constraint)
			constraints.append(And(r_constraints))
	return simplify(Or(constraints))
#exclude_graph for bound b returns constraints 
#telling solver to only find paths which all of their states
#(except for the initial and final state of the path) are not 
#included in the graph
def exclude_graph(graph, bound):

	#if the bound is one, only new edges are accepted
	if bound==1:
		constraints = []
		for e in graph.edge_list: 
			n1 = e.n1
			n2 = e.n2
			
			var_values = []
			for s in n1.var_dict: 
				x = Int(s + '.0')
				temp_const = (x == n1.var_dict[s])
				var_values.append(temp_const)

			for s in n2.var_dict: 
				x = Int(s + '.1')
				temp_const = (x == n2.var_dict[s])
				var_values.append(temp_const)

			constraints.append(And(var_values))
		return Not(Or(constraints))

	else: 
		constraints = []
		for n in graph.node_list: 
			for i in range(1, bound): 
				var_values = []
				for s in n.var_dict: 
					x = Int (s + '.' + str(i))
					temp_const = (x == n.var_dict[s])
					var_values.append(temp_const)
				constraints.append(And(var_values))
		return Not(Or(constraints))

def exclude_graph_2(graph, bound):

	#if the bound is one, only new edges are accepted
	if bound==1:
		constraints = []
		for e in graph.edge_list: 
			n1 = e.n1
			n2 = e.n2
			
			var_values = []
			for s in n1.var_dict: 
				x = Int(s + '.0')
				temp_const = (x == n1.var_dict[s])
				var_values.append(temp_const)

			for s in n2.var_dict: 
				x = Int(s + '.1')
				temp_const = (x == n2.var_dict[s])
				var_values.append(temp_const)

			constraints.append(And(var_values))
		return Not(Or(constraints))

	else: 
		constraints = []
		for n in graph.node_list: 
			for i in range(1, bound+1): 
				var_values = []
				for s in n.var_dict: 
					x = Int (s + '.' + str(i))
					temp_const = (x == n.var_dict[s])
					var_values.append(temp_const)
				constraints.append(And(var_values))
		return Not(Or(constraints))


#returns the set of contraints ensuring a path found of 
#specified bound does not contain a loop (two equal nodes does
#not appear in the path). Works incrementally: for ensuring a path
#of length 6 has no loops returned constraints should be added 
#to solver for bounds 4, 5, 6
#Assuming we are adding paths to the graph, and that we exclude the 
#graph beforehand in the solver, only paths with length >=4
#might have a loop
def loop_constraint(model, bound):
	if bound<4:
		return True
	else:
		species_vector = model.species_vector()
		constraints = []
		for s in species_vector:
			x = Int(s + '.' + str(bound-1))
			for i in range(1,bound-1): 
				y = Int(s + '.' + str(i))
				constraints.append(x==y)
		return Not(Or(constraints))

def construct_path(graph, model, model_name, prism, csl_prop, cp_bound, count_, prob_, prob_bound, property_var, property_val, mc_step):
	count = count_
	prob = prob_
	curr_bound = 1

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

	while curr_bound<cp_bound:
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
			flag = True
			path = solver.model()
			graph.add_path(path)
			solver.add(exclude_path(path))
			count = count + 1
			if (count%mc_step)==0: 
				prob = graph.model_check(model, model_name, prism, csl_prop)
				print('total number of paths so far: ' + str(count))
				print('probability= ' + str(prob))
				print('+++++')
			if prob>=prob_bound:
				return count, prob, True

		curr_bound = curr_bound + 1
		if flag: 
			curr_bound = 1

	return count, prob, False


def construct_path_3(graph, model, model_name, prism, csl_prop, cp_bound, count_, prob_, prob_bound, property_var, property_val, mc_step):
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
			graph.add_path(path)
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


def construct_path_2(graph, model, model_name, prism, csl_prop, cp_bound, count_, prob_, prob_bound, property_var, property_val, mc_step):
	count = count_
	prob = prob_
	curr_bound = 1

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

	while curr_bound<cp_bound:
		print(curr_bound)
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

		bound_count = 0
		while(solver.check(target_const)==sat):
			bound_count = bound_count + 1
			flag = True
			path = solver.model()
			graph.add_path(path)
			solver.add(exclude_path(path))
			count = count + 1
			if (count%mc_step)==0: 
				prob = graph.model_check(model, model_name, prism, csl_prop)
				print('total number of paths so far: ' + str(count))
				print('probability= ' + str(prob))
				print('+++++')
			if prob>=prob_bound:
				return count, prob, True
			if bound_count>1800:
				break

		curr_bound = curr_bound + 1
		# if flag: 
		# 	curr_bound = 1

	return count, prob, False





