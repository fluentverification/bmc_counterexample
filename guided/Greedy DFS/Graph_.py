from Utils import get_total_outgoing_rate
import copy
import subprocess


class Node:
    def __init__(self, var_values = None, is_target = False
                 , is_initial = False, outgoing_edges = None
                 , reachability_probability = 0, in_graph = False
                 , candidate_edges = None, parent = None, index = -1, bound = None):
        self.var_values = var_values
        self.is_target = is_target
        self.is_initial = is_initial
        self.outgoing_edges = outgoing_edges
        self.reachability_probability = reachability_probability
        self.in_graph = in_graph
        self.candidate_edges = candidate_edges
        self.parent = parent
        self.index = index
        self.bound = bound

    def get_var_values(self):
        return self.var_values

    def set_var_values(self, var_values):
        if self.var_values == None:
            self.var_values = []
        self.var_values = var_values

    def get_is_target(self):
        return self.is_target

    def set_is_target(self, is_target):
        self.is_target = is_target

    def get_is_initial(self):
        return self.is_initial

    def set_is_initial(self, is_initial):
        self.is_initial = is_initial

    def get_outgoing_edges(self):
        return self.outgoing_edges

    def add_outgoing_edge(self, outgoing_edge):
        if self.outgoing_edges == None:
            self.outgoing_edges = []
        self.outgoing_edges.append(outgoing_edge)

    def contain_outgoing_edge(self, edge):
        if self.get_outgoing_edges() == None:
            return False
        for e in self.get_outgoing_edges():
            if e.equals(edge):
                return True
        return False
    
    def get_reachability_probability(self):
        return self.reachability_probability

    def set_reachability_probability(self, reachability_probability):
        self.reachability_probability = reachability_probability

    def get_in_graph(self):
        return self.in_graph

    def set_in_graph(self, in_graph):
        self.in_graph = in_graph

    def get_candidate_edges(self):
        return self.candidate_edges

    def add_candidate_edge(self, candidate_edge):
        if self.candidate_edges == None:
            self.candidate_edges = []
        self.candidate_edges.append(candidate_edge)

    def equals(self, node):
        return self.var_values == node.get_var_values()
    
    def get_parent(self):
        return self.parent

    def set_parent(self, parent):
        self.parent = parent

    def get_index(self):
        return self.index
    
    def set_index(self, index):
        self.index = index

    def get_bound(self):
        return self.bound
    
    def set_bound(self, bound):
        self.bound = bound

class Edge:
    def __init__(self, src = None, dst = None
                 , in_graph = False, rate = 0, prob_dtmc = 0, reaction = -1):
        self.src = src
        self.dst = dst
        self.in_graph = in_graph
        self.rate = rate
        self.reaction = reaction
        self.prob_dtmc = prob_dtmc

    def get_src(self):
        return self.src
    
    def set_src(self, src):
        self.src = src

    def get_dst(self):
        return self.dst
    
    def set_dst(self, dst):
        self.dst = dst

    def get_in_graph(self):
        return self. in_graph
    
    def set_in_graph(self, in_graph):
        self.in_graph = in_graph

    def get_rate(self):
        return self.rate
    
    def set_rate(self, rate):
        self.rate = rate

    def get_prob_dtmc(self):
        return self.prob_dtmc
    
    def set_prob_dtmc(self, prob_dtmc):
        self.prob_dtmc = prob_dtmc

    def get_reaction(self):
        return self.reaction
    
    def set_reaction(self, reaction):
        self.reaction = reaction

    def equals(self, edge):
        if self.src.equals(edge.src):
            if self.dst.equals(edge.dst):
                if self.reaction == edge.reaction:
                    return True
        return False

class Graph:

    def __init__(self, node_list = None, edge_list = None, sink_state = None):
        self.node_list = node_list
        self.edge_list = edge_list
        self.sink_state = sink_state

    def get_node_list(self):
        return self.node_list
    
    def add_node(self, node):
        if self.node_list == None:
            self.node_list = []
        self.node_list.append(node)

    def get_edge_list(self):
        return self.edge_list
    
    def add_edge(self, edge):
        if self.edge_list == None:
            self.edge_list = []
        self.edge_list.append(edge)

    def get_size(self):
        if self.node_list != None:
            if self.edge_list != None:
                return len(self.node_list) + len(self.edge_list)
            else:
                return len(self.node_list)
        else:
            return 0

    def set_sink_state(self, node):
        self.add_node(node)
        self.sink_state = node

    def get_sink_state(self):
        return self.sink_state

    def contain_edge(self, edge):
        if self.get_edge_list() == None:
            return [False, None]
        for e in self.get_edge_list():
            if e.equals(edge):
                return [True, e]
        return [False, None]

    def contain_node(self, node):
        if self.get_node_list() == None:
            return [False, None]
        for n in self.get_node_list():
            if n.equals(node):
                return [True, n]
        return [False, None]
    
    def to_file_prism(self, file_name_prefix, model):
        # the states file (.sta)
        species_vector = model.species_vector()
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

            for i, n in enumerate(self.node_list): 
                if n.get_is_initial():
                #keep the initial state index for .lab file
                    initial_state_index = i

                line = str(i) + ':('
                for j in n.get_var_values(): 
                    line = line + str(j) + ','
                line_list = list(line)
                line_list[-1] = ')'
                line = ''.join(line_list)
                f.write(line)
                f.write('\n')
                n.set_index(i)

        # the labels file (.lab)
        labels_file_name = file_name_prefix + '.lab'
        with open(labels_file_name, mode='w', encoding='ascii') as f:
            f.truncate(0)
            lab_line = '0="init" 2="sink"'
            f.write(lab_line)
            f.write('\n')
            f.write(str(initial_state_index) + ': 0')
            f.close()

        # the transition file (.tra)
        trans_file_name = file_name_prefix + '.tra'
        with open(trans_file_name, mode='w', encoding='ascii') as f:
            f.truncate(0)
            size_line = str(len(self.node_list)) + ' ' + str(len(self.edge_list))
            f.write(size_line)
            f.write('\n')
            for e in self.edge_list: 
                line = str(e.get_src().get_index()) + ' ' + str(e.get_dst().get_index()) + ' ' + str(e.get_rate())
                f.write(line)
                f.write('\n')
            f.close()

    def add_sink_trans(self, model):
        graph = copy.deepcopy(self)
        sink_var_values = [-1] * len(model.species_vector())
        sink_node = Node(var_values=sink_var_values)

        for n in graph.get_node_list():
            current_outgoing_rate = 0
            if n.get_outgoing_edges() == None:
                continue
            for e in n.get_outgoing_edges():
                current_outgoing_rate = current_outgoing_rate + e.get_rate()
            total_outgoing_rate = get_total_outgoing_rate(n.get_var_values(), model)
            remaining_outgoing_rate = total_outgoing_rate - current_outgoing_rate
            sink_edge = Edge()
            sink_edge.set_src(n)
            sink_edge.set_dst(sink_node)
            sink_edge.set_rate(remaining_outgoing_rate)
            graph.add_edge(sink_edge)
        graph.add_node(sink_node)

        return graph
    
    def add_sink_trans_in_place(self, model):
        if self.get_sink_state() == None:
            sink_var_values = [-1] * len(model.species_vector())
            sink_node = Node(var_values=sink_var_values)
            self.set_sink_state(sink_node)

        sink_node = self.get_sink_state()

        for n in self.get_node_list():
            current_outgoing_rate = 0
            if n.get_outgoing_edges() == None:
                continue
            for e in n.get_outgoing_edges():
                current_outgoing_rate = current_outgoing_rate + e.get_rate()
            total_outgoing_rate = get_total_outgoing_rate(n.get_var_values(), model)
            remaining_outgoing_rate = total_outgoing_rate - current_outgoing_rate
            sink_edge = Edge()
            sink_edge.set_src(n)
            sink_edge.set_dst(sink_node)
            sink_edge.set_rate(remaining_outgoing_rate)
            self.add_edge(sink_edge)
    
    def model_check(self, model, prism_bin, csl_prop):
        new_graph = self.add_sink_trans(model)
        new_graph.to_file_prism("test", model)
        stdout_result = subprocess.run([prism_bin, '-importmodel', "test.all", '-pf', csl_prop, '-ctmc'], stdout=subprocess.PIPE)
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