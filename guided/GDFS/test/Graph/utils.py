from Utils.misc import get_total_outgoing_rate
import subprocess

def to_file_prism_format(graph, model, file_name_prefix):
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
            for n in graph.nodes.values():
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
            size_line = str(len(graph.nodes)) + ' ' + str(len(graph.edges))
            f.write(size_line)
            f.write('\n')
            for e in graph.edges.values(): 
                line = str(e.src.index) + ' ' + str(e.dst.index) + ' ' + str(e.rate)
                f.write(line)
                f.write('\n')

            for n in graph.nodes.values():
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

def check_probability(graph, model, file_name_prefix, prism_bin, csl_prop):
    to_file_prism_format(graph, model, "./results/" + file_name_prefix + "/" + file_name_prefix)
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