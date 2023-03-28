
node1 = Node()
node1.var_values = model.get_initial_state()
node1.initial_state = True

node2_var_values =  (1, 50, 0, 0, 49, 1)
node3_var_values =  (1, 51, 0, 1, 49, 0)
node4_var_values =  (1, 51, 0, 0, 48, 1)
node5_var_values =  (1, 52, 0, 1, 48, 0)
node6_var_values =  (1, 52, 0, 0, 47, 1)
node7_var_values =  (1, 53, 0, 1, 47, 0)
node8_var_values =  (1, 53, 0, 0, 46, 1)
node9_var_values =  (1, 54, 0, 1, 46, 0)
node10_var_values = (1, 54, 0, 0, 45, 1)
node11_var_values = (1, 55, 0, 1, 45, 0)
node12_var_values = (1, 55, 0, 0, 44, 1)
node13_var_values = (1, 56, 0, 1, 44, 0)
node14_var_values = (1, 56, 0, 0, 43, 1)
node15_var_values = (1, 57, 0, 1, 43, 0)
node16_var_values = (1, 57, 0, 0, 42, 1)
node17_var_values = (1, 58, 0, 1, 42, 0)
node18_var_values = (1, 58, 0, 0, 41, 1)
node19_var_values = (1, 59, 0, 1, 41, 0)
node20_var_values = (1, 59, 0, 0, 40, 1)

node2 = Node()
node2.var_values = node2_var_values
node3 = Node()
node3.var_values = node3_var_values
node4 = Node()
node4.var_values = node4_var_values
node5 = Node()
node5.var_values = node5_var_values
node6 = Node()
node6.var_values = node6_var_values
node7 = Node()
node7.var_values = node7_var_values
node8 = Node()
node8.var_values = node8_var_values
node9 = Node()
node9.var_values = node9_var_values
node10 = Node()
node10.var_values = node10_var_values
node11 = Node()
node11.var_values = node11_var_values
node12 = Node()
node12.var_values = node12_var_values
node13 = Node()
node13.var_values = node13_var_values
node14 = Node()
node14.var_values = node14_var_values
node15 = Node()
node15.var_values = node15_var_values
node16 = Node()
node16.var_values = node16_var_values
node17 = Node()
node17.var_values = node17_var_values
node18 = Node()
node18.var_values = node18_var_values
node19 = Node()
node19.var_values = node19_var_values
node20 = Node()
node20.var_values = node20_var_values

edge1 = Edge()
edge1.src = node1
edge1.dst = node2
edge1.reaction = 3
edge2 = Edge()
edge2.src = node2
edge2.dst = node3
edge2.reaction = 5
edge3 = Edge()
edge3.src = node3
edge3.dst = node4
edge3.reaction = 3
edge4 = Edge()
edge4.src = node4
edge4.dst = node5
edge4.reaction = 5
edge5 = Edge()
edge5.src = node5
edge5.dst = node6
edge5.reaction = 3
edge6 = Edge()
edge6.src = node6
edge6.dst = node7
edge6.reaction = 5
edge7 = Edge()
edge7.src = node7
edge7.dst = node8
edge7.reaction = 3
edge8 = Edge()
edge8.src = node8
edge8.dst = node9
edge8.reaction = 5
edge9 = Edge()
edge9.src = node9
edge9.dst = node10
edge9.reaction = 3
edge10 = Edge()
edge10.src = node10
edge10.dst = node11
edge10.reaction = 5
edge11 = Edge()
edge11.src = node11
edge11.dst = node12
edge11.reaction = 3
edge12 = Edge()
edge12.src = node12
edge12.dst = node13
edge12.reaction = 5
edge13 = Edge()
edge13.src = node13
edge13.dst = node14
edge13.reaction = 3
edge14 = Edge()
edge14.src = node14
edge14.dst = node15
edge14.reaction = 5
edge15 = Edge()
edge15.src = node15
edge15.dst = node16
edge15.reaction = 3
edge16 = Edge()
edge16.src = node16
edge16.dst = node17
edge16.reaction = 5
edge17 = Edge()
edge17.src = node17
edge17.dst = node18
edge17.reaction = 3
edge18 = Edge()
edge18.src = node18
edge18.dst = node19
edge18.reaction = 5
edge19 = Edge()
edge19.src = node19
edge19.dst = node20
edge19.reaction = 3

edge1.rate = get_reaction_rate(node1.var_values, model, reaction=3)
edge2.rate = get_reaction_rate(node2.var_values, model, reaction=5)
edge3.rate = get_reaction_rate(node3.var_values, model, reaction=3)
edge4.rate = get_reaction_rate(node4.var_values, model, reaction=5)
edge5.rate = get_reaction_rate(node5.var_values, model, reaction=3)
edge6.rate = get_reaction_rate(node6.var_values, model, reaction=5)
edge7.rate = get_reaction_rate(node7.var_values, model, reaction=3)
edge8.rate = get_reaction_rate(node8.var_values, model, reaction=5)
edge9.rate = get_reaction_rate(node9.var_values, model, reaction=3)
edge10.rate = get_reaction_rate(node10.var_values, model, reaction=5)
edge11.rate = get_reaction_rate(node11.var_values, model, reaction=3)
edge12.rate = get_reaction_rate(node12.var_values, model, reaction=5)
edge13.rate = get_reaction_rate(node13.var_values, model, reaction=3)
edge14.rate = get_reaction_rate(node14.var_values, model, reaction=5)
edge15.rate = get_reaction_rate(node15.var_values, model, reaction=3)
edge16.rate = get_reaction_rate(node16.var_values, model, reaction=5)
edge17.rate = get_reaction_rate(node17.var_values, model, reaction=3)
edge18.rate = get_reaction_rate(node18.var_values, model, reaction=5)
edge19.rate = get_reaction_rate(node19.var_values, model, reaction=3)



graph = Graph()

graph.add_edge(edge1)
graph.add_edge(edge2)
graph.add_edge(edge3)
graph.add_edge(edge4)
graph.add_edge(edge5)
graph.add_edge(edge6)
graph.add_edge(edge7)
graph.add_edge(edge8)
graph.add_edge(edge9)
graph.add_edge(edge10)
graph.add_edge(edge11)
graph.add_edge(edge12)
graph.add_edge(edge13)
graph.add_edge(edge14)
graph.add_edge(edge15)
graph.add_edge(edge16)
graph.add_edge(edge17)
graph.add_edge(edge18)
graph.add_edge(edge19)
