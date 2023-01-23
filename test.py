from z3 import *
import numpy as np
import models
from Graph import graph, node
from utils import *

constructor = getattr(models, 'circuit0x8E_010_111')
model = constructor()
X = model.initial_state()

V = model.reaction_rates(X)

for i in V:
	print(i)