from Utils import is_target, get_total_outgoing_rate, get_reaction_rate
import subprocess
import stormpy
import numpy as np
import scipy.sparse as sp


class Node:
    def __init__(self):
        self.var_values = tuple()
        self.initial_state = False
        self.out_edges = {}
        self.index = -1
        self.reachability_probability = 0
        self.in_edges = {}
        self.last_in_edge = None

    def add_out_edge(self, edge):
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.out_edges:
            self.out_edges[edge_tuple] = edge

    def add_in_edge(self, edge):
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.in_edges:
            self.in_edges[edge_tuple] = edge

    def __eq__(self, node) -> bool:
        return self.var_values == node.var_values

    def __lt__(self, node):
        return True
        
    

class Edge:
    def __init__(self):
        self.src = None
        self.dst = None
        self.rate = 0
        self.reaction = -1

    def __eq__(self, edge) -> bool:
        src_check = (self.src.var_values == edge.src.var_values)
        dst_check = (self.dst.var_values == edge.dst.var_values)
        reaction_check = (self.reaction == edge.reaction)
        return src_check and dst_check and reaction_check
    
    def get_tuple(self):
        return self.src.var_values + self.dst.var_values + (self.reaction,)



class Graph:
    def __init__(self):
        self.nodes = {}
        self.edges = {}

    def get_node(self, node_tuple):
        if node_tuple in self.nodes:
            return (True, self.nodes[node_tuple])
        return (False, None)
    
    def get_edge(self, edge_tuple):
        if edge_tuple in self.edges:
            return (True, self.edges[edge_tuple])
        return (False, None)
    
    def add_node(self, node):
        if node.var_values not in self.nodes:
            self.nodes[node.var_values] = node
            return node
        else:
            return self.nodes[node.var_values]

    def add_edge(self, edge):
        edge.src = self.add_node(edge.src)
        edge.dst = self.add_node(edge.dst)
        edge_tuple = edge.get_tuple()
        if edge_tuple not in self.edges:
            self.edges[edge_tuple] = edge
            self.nodes[edge.src.var_values].add_out_edge(edge)
            self.nodes[edge.dst.var_values].add_in_edge(edge)
    
    def to_file_prism_format(self, model, file_name_prefix):
        # the states file (.sta) ###########################################
        species_vector = model.get_species_tuple()
        states_file_name = file_name_prefix + '.sta'     
        with open(states_file_name, mode='w', encoding='ascii') as f:
            f.truncate(0)       
            #first line of the .sta file shows the order of variables
            state_vector_line = '('
            for e in species_vector: 
                state_vector_line = state_vector_line + e + ','
            state_vector_line_list = list(state_vector_line)
            state_vector_line_list[-1] = ')'
            state_vector_line = ''.join(state_vector_line_list)
            f.write(state_vector_line)
            f.write('\n')

            #every other line is a unique state
            index = 0
            for n in self.nodes.values():
                n.index = index
                if n.initial_state:
                    initial_state_index = index
                line = str(index) + ':('
                for j in n.var_values: 
                    line = line + str(j) + ','
                line_list = list(line)
                line_list[-1] = ')'
                line = ''.join(line_list)
                f.write(line)
                f.write('\n')
                index = index + 1
            #adding the sink state
            line = str(index) + ':('
            for i in range (len(model.get_species_tuple())):
                line = line + '-1,'
            line_list = list(line)
            line_list[-1] = ')'
            line = ''.join(line_list)
            f.write(line)
            f.write('\n')
            
        ##########################################################################

        # the labels file (.lab) #################################################
        labels_file_name = file_name_prefix + '.lab'
        with open(labels_file_name, mode='w', encoding='ascii') as f:
            f.truncate(0)
            lab_line = '0="init" 2="sink"'
            f.write(lab_line)
            f.write('\n')
            f.write(str(initial_state_index) + ': 0')
            f.write('\n')
            f.write(str(index) + ': 2')
            f.close()
        ##########################################################################

        # the transition file (.tra) #############################################
        trans_file_name = file_name_prefix + '.tra'
        with open(trans_file_name, mode='w', encoding='ascii') as f:
            f.truncate(0)
            size_line = str(len(self.nodes)) + ' ' + str(len(self.edges))
            f.write(size_line)
            f.write('\n')
            for e in self.edges.values(): 
                line = str(e.src.index) + ' ' + str(e.dst.index) + ' ' + str(e.rate)
                f.write(line)
                f.write('\n')

            for n in self.nodes.values():
                total_rate = get_total_outgoing_rate(n.var_values, model)
                current_rate = 0
                for e in n.out_edges.values():
                    current_rate = current_rate + e.rate
                remainin_rate = total_rate - current_rate
                if remainin_rate<=0:
                    continue
                line = str(n.index) + ' ' + str(index) + ' ' + str(remainin_rate)
                f.write(line)
                f.write('\n')
            f.close()

    def check_probability(self, model, file_name_prefix, prism_bin, csl_prop):
        self.to_file_prism_format(model, "./results/" + file_name_prefix + "/" + file_name_prefix)
        stdout_result = subprocess.run([prism_bin, '-importmodel', "./results/" + file_name_prefix + "/" + file_name_prefix + ".all", '-pf', csl_prop, '-ctmc'], stdout=subprocess.PIPE)
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

    #Uses STORM to model-check. Has errors!! debug before deploying
    def model_check(self, model, target_index, target_value, property):
        num_nodes = len(self.nodes)
        exit_rates = np.empty(num_nodes+1, dtype='float64')
        rate_matrix = sp.csr_matrix((num_nodes+1, num_nodes+1), dtype='float64')
        rate_matrix = rate_matrix.tolil()
        rate_matrix = np.zeros((num_nodes+1, num_nodes+1), dtype='float64')
        target_states = []
        init_states = []
        property = stormpy.parse_properties_without_context(property)[0]
    
        sink_node_tuple = tuple()
        for i in range(len(model.get_species_tuple())):
            sink_node_tuple = sink_node_tuple + (-1,)
        sink_node = Node()
        sink_node.var_values = sink_node_tuple
        sink_node.index = num_nodes
        self.add_node(sink_node)

        index = 0
        for n in self.nodes.values():
            n.index = index
            index = index + 1
            if n.initial_state:
                init_states.append(n.index)
            if is_target(n.var_values, target_index, target_value):
                target_states.append(n.index)
    
        if len(init_states) != 1:
            raise Exception('|number of initial states does not equal to 1|')
    
        for n in self.nodes.values():
            sum_out_rates = 0
            for e in n.out_edges.values():
                rate_matrix[n.index, e.dst.index] = e.rate
                sum_out_rates = sum_out_rates + e.rate
            total_outgoing_rate = get_total_outgoing_rate(n.var_values, model)
            exit_rates[n.index] = total_outgoing_rate
            remaining_rate = total_outgoing_rate - sum_out_rates
            if remaining_rate>0:
                rate_matrix[n.index, sink_node.index] = remaining_rate
    
        state_labeling = stormpy.storage.StateLabeling(num_nodes+1)
        state_labels = {'init', 'target'}
        for label in state_labels:
            state_labeling.add_label(label)

        for i in init_states:
            state_labeling.add_label_to_state('init', i)
    
        for i in target_states:
            state_labeling.add_label_to_state('target', i)
    
        rate_matrix = stormpy.build_sparse_matrix(rate_matrix)
        components = stormpy.SparseModelComponents(transition_matrix=rate_matrix, state_labeling=state_labeling, rate_transitions=True)
        components.exit_rates = exit_rates.tolist()
        ctmc = stormpy.storage.SparseCtmc(components)
        result = stormpy.model_checking(ctmc, property)
        return result.at(init_states[0])

