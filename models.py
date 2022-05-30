import numpy as np

#species cannot have dot "." in their name as it is used internally
#by the program to keep the state of the spcecies through time(steps)

class enzymatic_futile_cycle: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
		self.species_vector_ = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
		
		#same size as species_vector_
		self.initial_state_ = [1, 50, 0, 1, 50, 0]

		#each Ri_in and Ri_out is the same size as species_vector
		R1_in  = [1, 1, 0, 0, 0, 0]
		R1_out = [0, 0, 1, 0, 0, 0]

		R2_in  = [0, 0, 1, 0, 0, 0]
		R2_out = [1, 1, 0, 0, 0, 0]

		R3_in  = [0, 0, 1, 0, 0, 0]
		R3_out = [1, 0, 0, 0, 1, 0]

		R4_in  = [0, 0, 0, 1, 1, 0]
		R4_out = [0, 0, 0, 0, 0, 1]

		R5_in  = [0, 0, 0, 0, 0, 1]
		R5_out = [0, 0, 0, 1, 1, 0]

		R6_in  = [0, 0, 0, 0, 0, 1]
		R6_out = [0, 1, 0, 1, 0, 0]

		R = [[R1_in, R1_out], [R2_in, R2_out], [R3_in, R3_out], [R4_in, R4_out], [R5_in, R5_out], [R6_in, R6_out]]
		
		#same size as R
		self.reaction_rates_ = [1.0, 1.0, 0.1, 1.0, 1.0, 0.1]
		############# end of model descriptioon block #############
		###########################################################
		###########################################################

		#reactions is a dictionary. For each key i in {1, 2, 3, ...}
		#The value would be a vector of size 3: [Ri_in, Ri_out, Ri_effect]
		self.reactions = {}
		for i, e in enumerate(R):
			R_effect = np.array(np.array(e[1])-np.array(e[0])).tolist()
			self.reactions[i+1] = [e[0], e[1], R_effect] 
	
	def species_vector(self):
		return  self.species_vector_

	def initial_state(self):
		return self.initial_state_

	def reactions_vector(self): 
		return self.reactions

	def reaction_rates(self): 
		return self.reaction_rates_


		
###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################



class yeast_polarization: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
		self.species_vector_ = ['R', 'L', 'RL', 'G', 'Ga', 'Gbg', 'Gd']
		
		#same size as species_vector_
		self.initial_state_ = [50, 2, 0, 50, 0, 0, 0]

		#each Ri_in and Ri_out is the same size as species_vector
		#		 ['R', 'L', 'RL', 'G', 'Ga', 'Gbg', 'Gd']
		R1_in  = [0  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]
		R1_out = [1  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]

		R2_in  = [1  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]
		R2_out = [0  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]

		R3_in  = [1  ,  1 ,  0  ,  0 ,  0  ,  0   ,  0]
		R3_out = [0  ,  1 ,  1  ,  0 ,  0  ,  0   ,  0]

		R4_in  = [0  ,  0 ,  1  ,  0 ,  0  ,  0   ,  0]
		R4_out = [1  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]

		R5_in  = [0  ,  0 ,  1  ,  1 ,  0  ,  0   ,  0]
		R5_out = [0  ,  0 ,  0  ,  0 ,  1  ,  1   ,  0]

		R6_in  = [0  ,  0 ,  0  ,  0 ,  1  ,  0   ,  0]
		R6_out = [0  ,  0 ,  0  ,  0 ,  0  ,  0   ,  1]

		R7_in  = [0  ,  0 ,  0  ,  0 ,  0  ,  1   ,  1]
		R7_out = [0  ,  0 ,  0  ,  1 ,  0  ,  0   ,  0]

		R8_in  = [0  ,  0 ,  0  ,  0 ,  0  ,  0   ,  0]
		R8_out = [0  ,  0 ,  1  ,  0 ,  0  ,  0   ,  0]

		R = [[R1_in, R1_out], [R2_in, R2_out], [R3_in, R3_out], [R4_in, R4_out], [R5_in, R5_out], [R6_in, R6_out], [R7_in, R7_out], [R8_in, R8_out]]
		
		#same size as R
		self.reaction_rates_ = [0.0038, 4.0e-4, 0.042, 0.01, 0.011, 0.1, 1.05e3, 3.21]
		############# end of model descriptioon block #############
		###########################################################
		###########################################################

		#reactions is a dictionary. For each key i in {1, 2, 3, ...}
		#The value would be a vector of size 3: [Ri_in, Ri_out, Ri_effect]
		self.reactions = {}
		for i, e in enumerate(R):
			R_effect = np.array(np.array(e[1])-np.array(e[0])).tolist()
			self.reactions[i+1] = [e[0], e[1], R_effect] 
	
	def species_vector(self):
		return  self.species_vector_

	def initial_state(self):
		return self.initial_state_

	def reactions_vector(self): 
		return self.reactions

	def reaction_rates(self): 
		return self.reaction_rates_

###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################

class motility_regulation: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
		self.species_vector_ = ['codY', 'flache', 'SigD_hag', 'hag', 'CodY_flache', 'CodY_hag', 'CodY', 'SigD', 'Hag']
		
		#same size as species_vector_
		#                     ['codY', 'flache', 'SigD_hag', 'hag', 'CodY_flache', 'CodY_hag', 'CodY', 'SigD', 'Hag']
		self.initial_state_ = [  1   ,    1    ,      1    ,   1  ,       1      ,      1    ,  10   ,   10  ,   10 ]

		#each Ri_in and Ri_out is the same size as species_vector
		#          ['codY', 'flache', 'SigD_hag', 'hag', 'CodY_flache', 'CodY_hag', 'CodY', 'SigD', 'Hag']
		R1_in  =   [   1  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]
		R1_out  =  [   1  ,     0   ,     0     ,   0  ,      0       ,     0     ,   1   ,   0   ,   0  ]

		R2_in  =   [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   1   ,   0   ,   0  ]
		R2_out  =  [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]

		R3_in  =   [   0  ,     1   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]
		R3_out  =  [   0  ,     1   ,     0     ,   0  ,      0       ,     0     ,   0   ,   1   ,   0  ]

		R4_in  =   [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   1   ,   0  ]
		R4_out  =  [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]

		R5_in  =   [   0  ,     0   ,     1     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]
		R5_out  =  [   0  ,     0   ,     0     ,   1  ,      0       ,     0     ,   0   ,   1   ,   1  ]

		R6_in  =   [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   1  ]
		R6_out  =  [   0  ,     0   ,     0     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]

		R7_in  =   [   0  ,     0   ,     0     ,   1  ,      0       ,     0     ,   0   ,   1   ,   0  ]
		R7_out  =  [   0  ,     0   ,     1     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]

		R8_in  =   [   0  ,     0   ,     1     ,   0  ,      0       ,     0     ,   0   ,   0   ,   0  ]
		R8_out  =  [   0  ,     0   ,     0     ,   1  ,      0       ,     0     ,   0   ,   1   ,   0  ]

		R9_in  =   [   0  ,     1   ,     0     ,   0  ,      0       ,     0     ,   1   ,   0   ,   0  ]
		R9_out  =  [   0  ,     0   ,     0     ,   0  ,      1       ,     0     ,   0   ,   0   ,   0  ]

		R10_in  =  [   0  ,     0   ,     0     ,   0  ,      1       ,     0     ,   0   ,   0   ,   0  ]
		R10_out  = [   0  ,     1   ,     0     ,   0  ,      0       ,     0     ,   1   ,   0   ,   0  ]

		R11_in  =  [   0  ,     0   ,     0     ,   1  ,      0       ,     0     ,   1   ,   0   ,   0  ]
		R11_out  = [   0  ,     0   ,     0     ,   0  ,      0       ,     1     ,   0   ,   0   ,   0  ]

		R12_in  =  [   0  ,     0   ,     0     ,   0  ,      0       ,     1     ,   0   ,   0   ,   0  ]
		R12_out  = [   0  ,     0   ,     0     ,   1  ,      0       ,     0     ,   1   ,   0   ,   0  ]

		R = [[R1_in, R1_out], [R2_in, R2_out], [R3_in, R3_out], [R4_in, R4_out], [R5_in, R5_out], [R6_in, R6_out], [R7_in, R7_out], [R8_in, R8_out], [R9_in, R9_out], [R10_in, R10_out], [R11_in, R11_out], [R12_in, R12_out]]
		
		#same size as R
		#                       R1     R2     R3     R4    R5     R6    R7    R8    R9  R10  R11   R12
		self.reaction_rates_ = [0.1, 0.0002, 1.0, 0.0002, 1.0, 0.0002, 0.01, 0.1, 0.02, 0.1, 0.01, 0.1]
		############# end of model descriptioon block #############
		###########################################################
		###########################################################

		#reactions is a dictionary. For each key i in {1, 2, 3, ...}
		#The value would be a vector of size 3: [Ri_in, Ri_out, Ri_effect]
		self.reactions = {}
		for i, e in enumerate(R):
			R_effect = np.array(np.array(e[1])-np.array(e[0])).tolist()
			self.reactions[i+1] = [e[0], e[1], R_effect] 
	
	def species_vector(self):
		return  self.species_vector_

	def initial_state(self):
		return self.initial_state_

	def reactions_vector(self): 
		return self.reactions

	def reaction_rates(self): 
		return self.reaction_rates_
		
###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################

class single_species: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
		self.species_vector_ = ['S1', 'S2']
		
		#same size as species_vector_
		self.initial_state_ = [1, 40]

		#each Ri_in and Ri_out is the same size as species_vector
		R1_in  = [1, 0]
		R1_out = [1, 1]

		R2_in  = [0, 1]
		R2_out = [0, 0]

		

		R = [[R1_in, R1_out], [R2_in, R2_out]]
		
		#same size as R
		self.reaction_rates_ = [1.0, 0.025]
		############# end of model descriptioon block #############
		###########################################################
		###########################################################

		#reactions is a dictionary. For each key i in {1, 2, 3, ...}
		#The value would be a vector of size 3: [Ri_in, Ri_out, Ri_effect]
		self.reactions = {}
		for i, e in enumerate(R):
			R_effect = np.array(np.array(e[1])-np.array(e[0])).tolist()
			self.reactions[i+1] = [e[0], e[1], R_effect] 
	
	def species_vector(self):
		return  self.species_vector_

	def initial_state(self):
		return self.initial_state_

	def reactions_vector(self): 
		return self.reactions

	def reaction_rates(self): 
		return self.reaction_rates_