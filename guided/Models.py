import numpy as np

class test: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
		self.species_vector_ = ['S1', 'S2', 'S3']
		
		#same size as species_vector_
		self.initial_state_ = [1, 1, 1]

		#each Ri_in and Ri_out is the same size as species_vector
		R1_in  = [1, 0, 0]
		R1_out = [1, 1, 1]

		R2_in  = [1, 0, 0]
		R2_out = [1, 0, 1]

		

		R = [[R1_in, R1_out], [R2_in, R2_out]]
		
		#same size as R
		self.reaction_rates_ = [9.0, 1.0]
		############# end of model descriptioon block #############
		###########################################################
		###########################################################

		self.reactions = {}
		for i, e in enumerate(R):
			R_effect = np.array(np.array(e[1])-np.array(e[0])).tolist()
			self.reactions[i+1] = [e[0], e[1], R_effect] 
	
	def species_vector(self):
		return  self.species_vector_

	def initial_state(self):
		return self.initial_state_

	def reactions_dict(self): 
		return self.reactions

	def reaction_rates(self): 
		return self.reaction_rates_




###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################



#species cannot have dot "." in their name as it is used internally
#by the program to keep the state of the species through time(steps)

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

	def reactions_dict(self): 
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
		#R = [[R3_in, R3_out], [R5_in, R5_out], [R8_in, R8_out]]


		
		#same size as R
		self.reaction_rates_ = [0.0038, 4.0e-4, 0.042, 0.01, 0.011, 0.1, 1.05e3, 3.21]
		#self.reaction_rates_ = [0.042, 0.011, 3.21]

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

	def reactions_dict(self): 
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

	def reactions_dict(self): 
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




###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
###############################################################

class circuit0x8E: 
	def __init__(self): 
		
		###########################################################
		####### change this block for adding a new model ##########
																					
		self.species_vector_ = ['AmtR_protein', 									
		'Ara_AraC_protein', 
		'BetI_protein', 
		'HlyIIR_protein', 
		'LacI_protein', 
		'PhlF_protein', 
		'TetR_protein', 
		'YFP_protein',
		'topModel_AmtRpart_module_sub__pBAD', 
		'topModel_AmtRpart_module_sub__pHlyIIR', 
		'topModel_BetIpart_module_sub__pHlyIIR', 
		'topModel_BetIpart_module_sub__pTet',
		'topModel_HlyIIRpart_module_sub__pBAD', 
		'topModel_HlyIIRpart_module_sub__pTet', 
		'topModel_PhlFpart_module_sub__pAmtR', 
		'topModel_PhlFpart_module_sub__pTac', 
		'topModel_YFPpart_module_sub__pBetI', 
		'topModel_YFPpart_module_sub__pPhlF']

		#same size as species_vector_
		self.initial_state_ = [70,   # AmtR_protein         
      	60,   # Ara_AraC_protein     
      	70,   # BetI_protein         
      	0,    # HlyIIR_protein        
      	0,    # LacI_protein         
      	70,   # PhlF_protein         
      	0,    # TetR_protein        
      	0,    # YFP_protein          
      	2,    # topModel_AmtRpart_module_sub__pBAD      
      	2,    # topModel_AmtRpart_module_sub__pHlyIIR
      	2,    # topModel_BetIpart_module_sub__pHlyIIR
      	2,    # topModel_BetIpart_module_sub__pTet
      	2,    # topModel_HlyIIRpart_module_sub__pBAD
      	2,    # topModel_HlyIIRpart_module_sub__pTet
      	2,    # topModel_PhlFpart_module_sub__pAmtR
      	2,    # topModel_PhlFpart_module_sub__pTac
      	2,    # topModel_YFPpart_module_sub__pBetI
      	2]   # topModel_YFPpart_module_sub__pPhlF

		#each Ri_in and Ri_out is the same size as species_vector
		# R1  = AmtR_degradation_interaction
		# R2  = BetI_degradation_interaction
		# R3  = HlyIIR_degradation_interaction
		# R4  = PhlF_degradation_interaction
		# R5  = topModel_AmtRpart_module_sub__AmtR_protein_interaction_0
		# R6  = topModel_AmtRpart_module_sub__AmtR_protein_interaction_1
		# R7  = topModel_BetIpart_module_sub__BetI_protein_interaction_0
		# R8  = topModel_BetIpart_module_sub__BetI_protein_interaction_1
		# R9  = topModel_HlyIIRpart_module_sub__HlyIIR_protein_interaction_0
		# R10 = topModel_HlyIIRpart_module_sub__HlyIIR_protein_interaction_1
		# R11 = topModel_PhlFpart_module_sub__PhlF_protein_interaction_0
		# R12 = topModel_PhlFpart_module_sub__PhlF_protein_interaction_1
		# R13 = topModel_YFPpart_module_sub__YFP_protein_interaction_0
		# R14 = topModel_YFPpart_module_sub__YFP_protein_interaction_1
		# R15 = YFP_degradation_interaction

		R1_in   =   [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R1_out   =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R2_in   =   [0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R2_out   =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R3_in   =   [0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R3_out   =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R4_in   =   [0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R4_out   =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R5_in   =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R5_out   =  [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R6_in   =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R6_out   =  [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R7_in   =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R7_out   =  [0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R8_in   =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R8_out   =  [0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R9_in   =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R9_out   =  [0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R10_in  =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R10_out  =  [0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R11_in  =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R11_out  =  [0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R12_in  =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R12_out  =  [0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R13_in  =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R13_out  =  [0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R14_in  =   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R14_out  =  [0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

		R15_in  =   [0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		R15_out  =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


		R = [[R1_in, R1_out],
		[R2_in, R2_out],
		 [R3_in, R3_out],
		 [R4_in, R4_out],
		 [R5_in, R5_out],
		 [R6_in, R6_out],
		 [R7_in, R7_out],
		 [R8_in, R8_out],
		 [R9_in, R9_out],
		 [R10_in, R10_out],
		 [R11_in, R11_out],
		 [R12_in, R12_out],
		 [R13_in, R13_out],
		 [R14_in, R14_out],
		 [R15_in, R15_out]]
		
		
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

	def reactions_dict(self): 
		return self.reactions

	def reaction_rates(self):
		return [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
	
	def reaction_propensities(self, var_values): 
		X = var_values
				
		#Define ks:
		ka = 0.25       # Activated production rate
		ka_f = 0.0033   # Forward activation binding rate
		ka_r = 1.0      # Reverse activation binding rate
		kao_f = 1.0    # Forward activated RNAP binding rate
		kao_r = 1.0    # Reverse activated RNAP binding rate
		kb = 1.0E-4     # Basal production rate
		kc_f = 0.05    # Forward complex formation rate
		kc_r = 1.0     # Reverse complex formation rate
		kd = 7.5E-4    # Degradation rate 
		ko = 0.05      # Open complex production rate
		ko_f = 0.033    # Forward RNAP binding rate
		ko_r = 1.0     # Reverse RNAP binding rate
		kr_f = 0.5     # Forward repression binding rate
		kr_r = 1.0     # Reverse repression binding rate
		nc = 2.0       # Stoichiometry of binding
		ng = 2.0      # Initial promoter count
		np = 10.0       # Stoichiometry of production
		nr = 30.0      # Initial RNAP count

		topModel_AmtRpart_module_sub__ka = 0.25    
		topModel_AmtRpart_module_sub__ka_f = 0.0033 
		topModel_AmtRpart_module_sub__ka_r = 1
		topModel_AmtRpart_module_sub__kao_f = 1
		topModel_AmtRpart_module_sub__kao_r = 1
		topModel_AmtRpart_module_sub__kb = 1.0E-4
		topModel_AmtRpart_module_sub__ko = 0.05
		topModel_AmtRpart_module_sub__ko_f = 0.033
		topModel_AmtRpart_module_sub__ko_r = 1
		topModel_AmtRpart_module_sub__kr_f = 0.5
		topModel_AmtRpart_module_sub__kr_r = 1
		topModel_AmtRpart_module_sub__nc = 2
		topModel_AmtRpart_module_sub__ng = 2
		topModel_AmtRpart_module_sub__np = 10
		topModel_AmtRpart_module_sub__nr = 30
		topModel_BetIpart_module_sub__ka = 0.25
		topModel_BetIpart_module_sub__ka_f = 0.0033
		topModel_BetIpart_module_sub__ka_r = 1
		topModel_BetIpart_module_sub__kao_f = 1
		topModel_BetIpart_module_sub__kao_r = 1
		topModel_BetIpart_module_sub__kb = 1.0E-4
		topModel_BetIpart_module_sub__ko = 0.05
		topModel_BetIpart_module_sub__ko_f = 0.033
		topModel_BetIpart_module_sub__ko_r = 1
		topModel_BetIpart_module_sub__kr_f = 0.5
		topModel_BetIpart_module_sub__kr_r = 1
		topModel_BetIpart_module_sub__nc = 2;
		topModel_BetIpart_module_sub__ng = 2
		topModel_BetIpart_module_sub__np = 10
		topModel_BetIpart_module_sub__nr = 30
		topModel_HlyIIRpart_module_sub__ka = 0.25
		topModel_HlyIIRpart_module_sub__ka_f = 0.0033
		topModel_HlyIIRpart_module_sub__ka_r = 1
		topModel_HlyIIRpart_module_sub__kao_f = 1
		topModel_HlyIIRpart_module_sub__kao_r = 1
		topModel_HlyIIRpart_module_sub__kb = 1.0E-4
		topModel_HlyIIRpart_module_sub__ko = 0.05
		topModel_HlyIIRpart_module_sub__ko_f = 0.033
		topModel_HlyIIRpart_module_sub__ko_r = 1
		topModel_HlyIIRpart_module_sub__kr_f = 0.5
		topModel_HlyIIRpart_module_sub__kr_r = 1
		topModel_HlyIIRpart_module_sub__nc = 2
		topModel_HlyIIRpart_module_sub__ng = 2
		topModel_HlyIIRpart_module_sub__np = 10
		topModel_HlyIIRpart_module_sub__nr = 30
		topModel_PhlFpart_module_sub__ka = 0.25
		topModel_PhlFpart_module_sub__ka_f = 0.0033
		topModel_PhlFpart_module_sub__ka_r = 1
		topModel_PhlFpart_module_sub__kao_f = 1
		topModel_PhlFpart_module_sub__kao_r = 1
		topModel_PhlFpart_module_sub__kb = 1.0E-4
		topModel_PhlFpart_module_sub__ko = 0.05
		topModel_PhlFpart_module_sub__ko_f = 0.033
		topModel_PhlFpart_module_sub__ko_r = 1
		topModel_PhlFpart_module_sub__kr_f = 0.5
		topModel_PhlFpart_module_sub__kr_r = 1
		topModel_PhlFpart_module_sub__nc = 2
		topModel_PhlFpart_module_sub__ng = 2
		topModel_PhlFpart_module_sub__np = 10
		topModel_PhlFpart_module_sub__nr = 30
		topModel_YFPpart_module_sub__ka = 0.25
		topModel_YFPpart_module_sub__ka_f = 0.0033
		topModel_YFPpart_module_sub__ka_r = 1
		topModel_YFPpart_module_sub__kao_f = 1
		topModel_YFPpart_module_sub__kao_r = 1
		topModel_YFPpart_module_sub__kb = 1.0E-4
		topModel_YFPpart_module_sub__ko = 0.05
		topModel_YFPpart_module_sub__ko_f = 0.033
		topModel_YFPpart_module_sub__ko_r = 1
		topModel_YFPpart_module_sub__kr_f = 0.5
		topModel_YFPpart_module_sub__kr_r = 1
		topModel_YFPpart_module_sub__nc = 2
		topModel_YFPpart_module_sub__ng = 2
		topModel_YFPpart_module_sub__np = 10
		topModel_YFPpart_module_sub__nr = 30

		# Reaction Propensity constants
		R1_r  = kd*X[1-1] # AmtR_degradation_interaction
		R2_r  = kd*X[3-1] # BetI_degradation_interaction
		R3_r  = kd*X[4-1] # HlyIIR_degradation_interaction
		R4_r  = kd*X[6-1] # PhlF_degradation_interaction
		R5_r  = X[9-1]*(topModel_AmtRpart_module_sub__kb*topModel_AmtRpart_module_sub__ko_f/topModel_AmtRpart_module_sub__ko_r*topModel_AmtRpart_module_sub__nr+topModel_AmtRpart_module_sub__ka*topModel_AmtRpart_module_sub__kao_f/topModel_AmtRpart_module_sub__kao_r*topModel_AmtRpart_module_sub__nr*(topModel_AmtRpart_module_sub__ka_f/topModel_AmtRpart_module_sub__ka_r*X[2-1])**topModel_AmtRpart_module_sub__nc)/(1+topModel_AmtRpart_module_sub__ko_f/topModel_AmtRpart_module_sub__ko_r*topModel_AmtRpart_module_sub__nr+topModel_AmtRpart_module_sub__kao_f/topModel_AmtRpart_module_sub__kao_r*topModel_AmtRpart_module_sub__nr*(topModel_AmtRpart_module_sub__ka_f/topModel_AmtRpart_module_sub__ka_r*X[2-1])**topModel_AmtRpart_module_sub__nc) # topModel_AmtRpart_module_sub__X(1)_interaction_0
		R6_r  = X[10-1]*topModel_AmtRpart_module_sub__ko*topModel_AmtRpart_module_sub__ko_f/topModel_AmtRpart_module_sub__ko_r*topModel_AmtRpart_module_sub__nr/(1+topModel_AmtRpart_module_sub__ko_f/topModel_AmtRpart_module_sub__ko_r*topModel_AmtRpart_module_sub__nr+(topModel_AmtRpart_module_sub__kr_f/topModel_AmtRpart_module_sub__kr_r*X[4-1])**topModel_AmtRpart_module_sub__nc) # topModel_AmtRpart_module_sub__X(1)_interaction_1
		R7_r  = X[11-1]*topModel_BetIpart_module_sub__ko*topModel_BetIpart_module_sub__ko_f/topModel_BetIpart_module_sub__ko_r*topModel_BetIpart_module_sub__nr/(1+topModel_BetIpart_module_sub__ko_f/topModel_BetIpart_module_sub__ko_r*topModel_BetIpart_module_sub__nr+(topModel_BetIpart_module_sub__kr_f/topModel_BetIpart_module_sub__kr_r*X[4-1])**topModel_BetIpart_module_sub__nc) # topModel_BetIpart_module_sub__X(3)_interaction_0
		R8_r  = X[12-1]*topModel_BetIpart_module_sub__ko*topModel_BetIpart_module_sub__ko_f/topModel_BetIpart_module_sub__ko_r*topModel_BetIpart_module_sub__nr/(1+topModel_BetIpart_module_sub__ko_f/topModel_BetIpart_module_sub__ko_r*topModel_BetIpart_module_sub__nr+(topModel_BetIpart_module_sub__kr_f/topModel_BetIpart_module_sub__kr_r*X[7-1])**topModel_BetIpart_module_sub__nc) # topModel_BetIpart_module_sub__X(3)_interaction_1
		R9_r  = X[14-1]*(topModel_HlyIIRpart_module_sub__kb*topModel_HlyIIRpart_module_sub__ko_f/topModel_HlyIIRpart_module_sub__ko_r*topModel_HlyIIRpart_module_sub__nr+topModel_HlyIIRpart_module_sub__ka*topModel_HlyIIRpart_module_sub__kao_f/topModel_HlyIIRpart_module_sub__kao_r*topModel_HlyIIRpart_module_sub__nr*(topModel_HlyIIRpart_module_sub__ka_f/topModel_HlyIIRpart_module_sub__ka_r*X[2-1])**topModel_HlyIIRpart_module_sub__nc)/(1+topModel_HlyIIRpart_module_sub__ko_f/topModel_HlyIIRpart_module_sub__ko_r*topModel_HlyIIRpart_module_sub__nr+topModel_HlyIIRpart_module_sub__kao_f/topModel_HlyIIRpart_module_sub__kao_r*topModel_HlyIIRpart_module_sub__nr*(topModel_HlyIIRpart_module_sub__ka_f/topModel_HlyIIRpart_module_sub__ka_r*X[2-1])**topModel_HlyIIRpart_module_sub__nc) # topModel_HlyIIRpart_module_sub__HlyIIR_protein_interaction_0
		R10_r = X[11-1]*topModel_HlyIIRpart_module_sub__ko*topModel_HlyIIRpart_module_sub__ko_f/topModel_HlyIIRpart_module_sub__ko_r*topModel_HlyIIRpart_module_sub__nr/(1+topModel_HlyIIRpart_module_sub__ko_f/topModel_HlyIIRpart_module_sub__ko_r*topModel_HlyIIRpart_module_sub__nr+(topModel_HlyIIRpart_module_sub__kr_f/topModel_HlyIIRpart_module_sub__kr_r*X[7-1])**topModel_HlyIIRpart_module_sub__nc) # topModel_HlyIIRpart_module_sub__HlyIIR_protein_interaction_1
		R11_r = X[16-1]*ko*ko_f/ko_r*nr/(1+ko_f/ko_r*nr+(kr_f/kr_r*X[5-1])**nc) # topModel_PhlFpart_module_sub__X(6)_interaction_0                                                                                                     
		R12_r = X[15-1]*topModel_PhlFpart_module_sub__ko*topModel_PhlFpart_module_sub__ko_f/topModel_PhlFpart_module_sub__ko_r*topModel_PhlFpart_module_sub__nr/(1+topModel_PhlFpart_module_sub__ko_f/topModel_PhlFpart_module_sub__ko_r*topModel_PhlFpart_module_sub__nr+(topModel_PhlFpart_module_sub__kr_f/topModel_PhlFpart_module_sub__kr_r*X[1-1])**topModel_PhlFpart_module_sub__nc) # topModel_PhlFpart_module_sub__X(6)_interaction_1                                                                                                          
		R13_r = X[18-1]*topModel_YFPpart_module_sub__ko*topModel_YFPpart_module_sub__ko_f/topModel_YFPpart_module_sub__ko_r*topModel_YFPpart_module_sub__nr/(1+topModel_YFPpart_module_sub__ko_f/topModel_YFPpart_module_sub__ko_r*topModel_YFPpart_module_sub__nr+(topModel_YFPpart_module_sub__kr_f/topModel_YFPpart_module_sub__kr_r*X[6-1])**topModel_YFPpart_module_sub__nc) # topModel_YFPpart_module_sub__YFP_protein_interaction_0                                                                                                        
		R14_r = X[17-1]*topModel_YFPpart_module_sub__ko*topModel_YFPpart_module_sub__ko_f/topModel_YFPpart_module_sub__ko_r*topModel_YFPpart_module_sub__nr/(1+topModel_YFPpart_module_sub__ko_f/topModel_YFPpart_module_sub__ko_r*topModel_YFPpart_module_sub__nr+(topModel_YFPpart_module_sub__kr_f/topModel_YFPpart_module_sub__kr_r*X[3-1])**topModel_YFPpart_module_sub__nc) # topModel_YFPpart_module_sub__YFP_protein_interaction_1                                                                                                         
		R15_r = kd*X[8-1] # YFP_degradation_interaction                                                                                                        
		

		#                       		R1    R2    R3    R4    R5    R6    R7    R8    R9    R10    R11    R12    R13    R14    R15
		self.reaction_propensities_ = [R1_r, R2_r, R3_r, R4_r, R5_r, R6_r, R7_r, R8_r, R9_r, R10_r, R11_r, R12_r, R13_r, R14_r, R15_r]
		return self.reaction_propensities_