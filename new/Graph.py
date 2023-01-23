import json
import numpy as np
import math
import subprocess

class node:
	def __init__(self, var_dict): 
		#a dictionary assining values to all of model's variables
		self.var_dict = var_dict	
		self.is_initial = False
		#index is used for keeping .sta and .tra files synched
		self.index = -1
		#a list containing all the edges going out of the node
		self.out_edges = []

	def make_initial(self): 
		self.is_initial = True

	def equals(self, node_): 
		if self.var_dict == node_.var_dict: 
			return True
		return False

	def set_index(self, index_): 
		self.index = index_

	def add_out_edge(self, edge_): 
		self.out_edges.append(edge_)


#an edge has two nodes n1:source node, n2: destination node
class edge: 
	def __init__(self, n1, n2, reaction): 
		self.n1 = n1
		self.n2 = n2
		#a reaction index in {1, 2, 3, ...}
		self.reaction = reaction

	def equals(self, edge_): 
		if self.n1.equals(edge_.n1):
			if self.n2.equals(edge_.n2): 
				if self.reaction == edge_.reaction:
					return True
		return False

	def get_node_list(self): 
		return [self.n1, self.n2]


class graph: 
	def __init__(self):
		self.node_list = []
		self.edge_list = []

	def add_node(self, node_): 
		flag = True
		for j in self.node_list:
			if node_.equals(j): 
				flag = False
				return_node = j
		if flag:
			self.node_list.append(node_)
			return_node = node_
		return return_node

	#if the edge_ is not already in the graph, adds the edge to the graph
	#updates the source node's list of outgoing edges
	def add_edge(self, edge_): 
		flag = True
		for j in self.edge_list: 
			if edge_.equals(j): 
				flag = False
		if flag:
			contains_node1 = False
			contains_node2 = False
			for n in self.node_list: 
				if n.equals(edge_.get_node_list()[0]): 
					contains_node1=True
					edge_.n1 = n
				if n.equals(edge_.get_node_list()[1]):
					contains_node2=True
					edge_.n2 = n
			
			if not (contains_node1): 
				self.node_list.append(edge_.get_node_list()[0])

			if not (contains_node2): 
				self.node_list.append(edge_.get_node_list()[1]) 
			
			self.edge_list.append(edge_)
			node1 = edge_.get_node_list()[0]
			node1.add_out_edge(edge_)


	def add_path(self, path, reaction_dict): 
		#find the last index (length) of the path
		index = -1
		for i in path.decls(): 
			index_temp = int(i.name()[i.name().rfind('.')+1:])
			if index_temp > index: 
				index = index_temp

		#a path is an ordered list of states(var assignments)
		node_ordered_list = [None] * (index+1)

		for i in path.decls():
			if 'rate' in i.name(): 
				continue
			index_temp = int(i.name()[i.name().rfind('.')+1:])
			if node_ordered_list[index_temp] == None: 
				node_ordered_list[index_temp] = []
			node_ordered_list[index_temp].append([i.name()[:i.name().find('.')], int(str(path[i]))])

		for i in range(index): 
		
			var_dict_from = {}
			for e in node_ordered_list[i]: 
				if e[0] != 'selected_reaction':
					var_dict_from[e[0]]=e[1]
				else: 
					selected_reaction = e[1]
			
			var_dict_to = {}
			for e in node_ordered_list[i+1]: 
				if e[0] != 'selected_reaction':
					var_dict_to[e[0]]=e[1]

			from_node = node(var_dict_from)
			to_node = node(var_dict_to)
			if selected_reaction == 1:
				reaction_dict["R1"] = reaction_dict["R1"] + 1
			elif selected_reaction == 2:
				reaction_dict["R2"] = reaction_dict["R2"] + 1
			elif selected_reaction == 3:
				reaction_dict["R3"] = reaction_dict["R3"] + 1
			elif selected_reaction == 4:
				reaction_dict["R4"] = reaction_dict["R4"] + 1
			elif selected_reaction == 5:
				reaction_dict["R5"] = reaction_dict["R5"] + 1
			elif selected_reaction == 6:
				reaction_dict["R6"] = reaction_dict["R6"] + 1
			
			edge_ = edge(from_node, to_node, selected_reaction)
			self.add_edge(edge_)


	
	def add_path_o(self, path): 
		#find the last index (length) of the path
		index = -1
		for i in path.decls(): 
			index_temp = int(i.name()[i.name().rfind('.')+1:])
			if index_temp > index: 
				index = index_temp

		#a path is an ordered list of states(var assignments)
		node_ordered_list = [None] * (index+1)

		for i in path.decls():
			if 'rate' in i.name(): 
				continue
			index_temp = int(i.name()[i.name().rfind('.')+1:])
			if node_ordered_list[index_temp] == None: 
				node_ordered_list[index_temp] = []
			node_ordered_list[index_temp].append([i.name()[:i.name().find('.')], int(str(path[i]))])

		for i in range(index): 
		
			var_dict_from = {}
			for e in node_ordered_list[i]: 
				if e[0] != 'selected_reaction':
					var_dict_from[e[0]]=e[1]
				else: 
					selected_reaction = e[1]
			
			var_dict_to = {}
			for e in node_ordered_list[i+1]: 
				if e[0] != 'selected_reaction':
					var_dict_to[e[0]]=e[1]

			from_node = node(var_dict_from)
			to_node = node(var_dict_to)
			edge_ = edge(from_node, to_node, selected_reaction)
			self.add_edge(edge_)
	#save the graph in PRISM readable format
	def to_file(self, file_name_prefix, model): 
		species_vector = model.species_vector()
		initial_state_index = -1
		sink_state_index = -2

		#the states file (.sta)
		states_file_name = file_name_prefix + '.sta'
		with open(states_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			
			#first line of the .sta file shows the variable formats
			state_vector_line = '('
			for e in species_vector: 
				state_vector_line = state_vector_line + e + ','
			state_vector_line_list = list(state_vector_line)
			state_vector_line_list[-1] = ')'
			state_vector_line = ''.join(state_vector_line_list)
			f.write(state_vector_line)
			f.write('\n')

			for i, n in enumerate(self.node_list): 
				if n.is_initial:
					#keep the initial state index for .lab file
					initial_state_index = i
				#set an index for each node for synch between .sta and .tra
				n.set_index(i)
				line = str(i) + ':('
				for s in species_vector: 
					line = line + str(n.var_dict[s]) + ','
				line_list = list(line)
				line_list[-1] = ')'
				line = ''.join(line_list)
				f.write(line)
				f.write('\n')
			#sink state has the largest index. It also has the value
			# (-1, -1, -1, ...) assigned to its variables
			sink_line = str(self.node_list[-1].index + 1) + ':('
			for e in species_vector:
				sink_line = sink_line + '-1,'
			sink_line_list = list(sink_line)
			sink_line_list[-1] = ')'
			sink_line = ''.join(sink_line_list)
			f.write(sink_line)
			f.close()

		#the labels file (.lab)
		labels_file_name = file_name_prefix + '.lab'
		with open(labels_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			sink_state_index = self.node_list[-1].index + 1
			lab_line = '0="init" 2="sink"'
			f.write(lab_line)
			f.write('\n')
			f.write(str(initial_state_index) + ': 0')
			f.write('\n')
			f.write(str(sink_state_index) + ': 2')
			f.close()

		#the transitions file (.tra)
		
		rate = float()
		#list of all the transitions already in the graph + added sink transitions
		#each trans_list entry is a vector of size three: [src, dest, rate]
		trans_list  = [[]] * len(self.node_list)
			
		for i, n in enumerate(self.node_list):
			#get the state-vector for node n (to compute rate based on that)
			var_values = [None] * len(species_vector)
			for j, s in enumerate(species_vector):
				var_values[j] = n.var_dict[s]

			#sum of the rate of all the outgoing transitions out of the node in the current graph
			current_rate = 0
			#sum of the rate of all the outgoing transitions out of the state in the original model
			total_rate = 0

			trans_list_temp = []
			for e in n.out_edges:
				comb = 1
				
				# 2R1 + R2 --K--> R3. the rate would be: #(R1) * #(R1) * #(R2) * K
				for it, r in enumerate(model.reactions_vector()[e.reaction][0]):
					for c in range(r): 
						comb = comb * var_values[it]
				rate = model.reaction_rates()[e.reaction-1]
				rate = rate * comb
				
				trans_list_temp.append([n.index, e.get_node_list()[1].index, rate])
				current_rate = current_rate + rate

			#all the reactions could possibly go out of the state(except for zero reactants (handled by comb*0))
			for j, r in enumerate(model.reactions_vector().values()):
				comb = 1
				rate = model.reaction_rates()[j]
				for it, v in enumerate(r[0]):
					for c in range(v): 
						comb = comb * var_values[it]
				rate = rate * comb
				total_rate = total_rate + rate

			remaining_rate = total_rate - current_rate
			if remaining_rate > 0:
				trans_list_temp.append([n.index, sink_state_index, remaining_rate])
			trans_list[i] = trans_list_temp

		trans_file_name = file_name_prefix + '.tra'
		with open(trans_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			num_transitions = 0
			for n in trans_list: 
				for e in n: 
					num_transitions = num_transitions + 1
			
			#number of nodes is the size of node_list + 1(sink state)
			first_line = str(len(self.node_list)+1) + ' ' + str(num_transitions)
			f.write(first_line)
			f.write('\n')
			for n in trans_list: 
				for e in n:
					line = str(int(e[0])) + ' ' + str(int(e[1])) + ' ' + str(e[2])
					f.write(line)
					f.write('\n')
			f.close()


	#save the graph in PRISM readable format
	def to_file_2(self, file_name_prefix, model): 
		species_vector = model.species_vector()
		initial_state_index = -1
		sink_state_index = -2

		#the states file (.sta)
		states_file_name = file_name_prefix + '.sta'
		with open(states_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			
			#first line of the .sta file shows the variable formats
			state_vector_line = '('
			for e in species_vector: 
				state_vector_line = state_vector_line + e + ','
			state_vector_line_list = list(state_vector_line)
			state_vector_line_list[-1] = ')'
			state_vector_line = ''.join(state_vector_line_list)
			f.write(state_vector_line)
			f.write('\n')

			for i, n in enumerate(self.node_list): 
				if n.is_initial:
					#keep the initial state index for .lab file
					initial_state_index = i
				#set an index for each node for synch between .sta and .tra
				n.set_index(i)
				line = str(i) + ':('
				for s in species_vector: 
					line = line + str(n.var_dict[s]) + ','
				line_list = list(line)
				line_list[-1] = ')'
				line = ''.join(line_list)
				f.write(line)
				f.write('\n')
			#sink state has the largest index. It also has the value
			# (-1, -1, -1, ...) assigned to its variables
			sink_line = str(self.node_list[-1].index + 1) + ':('
			for e in species_vector:
				sink_line = sink_line + '-1,'
			sink_line_list = list(sink_line)
			sink_line_list[-1] = ')'
			sink_line = ''.join(sink_line_list)
			f.write(sink_line)
			f.close()

		#the labels file (.lab)
		labels_file_name = file_name_prefix + '.lab'
		with open(labels_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			sink_state_index = self.node_list[-1].index + 1
			lab_line = '0="init" 2="sink"'
			f.write(lab_line)
			f.write('\n')
			f.write(str(initial_state_index) + ': 0')
			f.write('\n')
			f.write(str(sink_state_index) + ': 2')
			f.close()

		#the transitions file (.tra)
		
		rate = float()
		#list of all the transitions already in the graph + added sink transitions
		#each trans_list entry is a vector of size three: [src, dest, rate]
		trans_list  = [[]] * len(self.node_list)
			
		for i, n in enumerate(self.node_list):
			#get the state-vector for node n (to compute rate based on that)
			var_values = [None] * len(species_vector)
			for j, s in enumerate(species_vector):
				var_values[j] = n.var_dict[s]

			#sum of the rate of all the outgoing transitions out of the node in the current graph
			current_rate = 0
			#sum of the rate of all the outgoing transitions out of the state in the original model
			total_rate = 0

			trans_list_temp = []
			for e in n.out_edges:
				comb = 1
				
				# 2R1 + R2 --K--> R3. the rate would be: #(R1) * #(R1) * #(R2) * K
				# for it, r in enumerate(model.reactions_vector()[e.reaction][0]):
				# 	for c in range(r): 
				# 		comb = comb * var_values[it]
				rate = model.reaction_rates(var_values)[e.reaction-1]
				#rate = rate * comb
				
				trans_list_temp.append([n.index, e.get_node_list()[1].index, rate])
				current_rate = current_rate + rate

			#all the reactions could possibly go out of the state(except for zero reactants (handled by comb*0))
			for j, r in enumerate(model.reactions_vector().values()):
				comb = 1
				rate = model.reaction_rates(var_values)[j]
				# for it, v in enumerate(r[0]):
				# 	for c in range(v): 
				# 		comb = comb * var_values[it]
				rate = rate * comb
				total_rate = total_rate + rate

			remaining_rate = total_rate - current_rate
			if remaining_rate > 0:
				trans_list_temp.append([n.index, sink_state_index, remaining_rate])
			trans_list[i] = trans_list_temp

		trans_file_name = file_name_prefix + '.tra'
		with open(trans_file_name, mode='w', encoding='ascii') as f:
			f.truncate()
			num_transitions = 0
			for n in trans_list: 
				for e in n: 
					num_transitions = num_transitions + 1
			
			#number of nodes is the size of node_list + 1(sink state)
			first_line = str(len(self.node_list)+1) + ' ' + str(num_transitions)
			f.write(first_line)
			f.write('\n')
			for n in trans_list: 
				for e in n:
					line = str(int(e[0])) + ' ' + str(int(e[1])) + ' ' + str(e[2])
					f.write(line)
					f.write('\n')
			f.close()

	#calculate the probability of the graph with respect to
	#csl_prop
	def model_check(self, model, model_name, prism_bin, csl_prop):
		self.to_file('./results/' + model_name + '/' + model_name, model)
		prism_model_files = 'results/' + model_name + '/' + model_name + '.all'
		stdout_result = subprocess.run([prism_bin, '-importmodel', prism_model_files, '-pf', csl_prop, '-ctmc'], stdout=subprocess.PIPE)
		stdout_result = stdout_result.stdout.decode('utf-8')
		stdout_result = stdout_result.splitlines()

		result = ''
		for r in stdout_result:
			if 'Result' in r: 
				result = r
		result = result[result.rfind(':')+2:]
		if ' ' in result: 
			result = result[:result.find(' ')]
		result = float(result)
		return result
