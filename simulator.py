## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## notes:
##	1. we should simulate tissue specificity, as the modeling takes consideration of that;
##	2. when we say a factor, it's latent factor condensed from some original variables;
##	3. TODO: in next stage of simulation, we can use the true genotype with fake beta to generate the expression profile, as the genotype distribution (MAF) may not be exactly what we assume here;
##	4. xxx


import numpy as np


##=====================
##==== global variables
## notes: TODO
##	1. shall we simulate Spike and Slab?
##	2. xxx
##=====================
#==== chromosome
L = 0								# unit: basepair

#==== individual
n_individual = 0

#==== tissue
n_tissue = 0

#==== SNP
n_SNP = 0
SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
SNP_pos_list = []
SNP_beta_rep = {} #{tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess

#==== gene
n_gene = 0
gene_rep = {} #{individual:{tissue:[], ...}, ...}		# gene expression list for tissue samples for individuals
gene_pos_list = []

#==== SNP gene pos map
pos_map = {}							# {gene:[snp1, snp2], ...}

#==== factor_cell
n_factor_cell = 0
factor_cell_beta_rep = {} #{cell factor:[], ...}		# coefficient lists for cell factors
beta_factor_cell_rep = {} #{tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues

#==== factor_batch (analogous to cell factor pathway)
n_batch = 0							# n_batch = n_batch_individual + n_batch_sample
n_batch_individual = 0
n_batch_sample = 0
n_factor_batch = 0						# the number of latent batch factors
batch_individual_rep = {} #{individual:[], ...}			# batch variable lists for individuals
batch_tissue_sample_rep = {} #{individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
factor_batch_beta_rep = {} #{batch factor:[], ...}		# coefficient lists for batch factors
beta_factor_batch_rep = {} #{gene:[], ...}			# coefficient lists of batch factors for genes




##=========================
##==== simulating variables
## notes: TODO
##	1. simulate variables observed and latent excluding those along the generation process
##	2. xxx
##=========================
def simu_genotype():	# variables
	# to fill in:
	##SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	##SNP_pos_list = []

	global SNP_rep
	global SNP_pos_list
	global L

	##==== SNP_rep (uniform randomly draw MAF from [0,1))
	for individual in SNP_rep:
		for i in range(len(SNP_rep[individual])):
			# random draw from [0, 1)
			dosage = np.random.random_sample()
			SNP_rep[individual][i] = dosage

	##==== SNP_pos_list (uniform randomly draw POS from [0,1) * L)
	temp_rep = {}
	temp_list = []
	for i in range(len(SNP_pos_list)):
		while 1:
			pos = int(np.random.random_sample() * L)
			if pos in temp_rep:
				continue
			else:
				temp_rep[pos] = 1
				temp_list.append(pos)
				break
	temp_list = np.array(temp_list)
	temp_list = np.sort(temp_list)
	SNP_pos_list = temp_list
	return


def simu_gene():
	# to fill in:
	##gene_pos_list = []

	global gene_pos_list

	##==== gene_pos_list
	temp_rep = {}
	temp_list = []
	for i in range(len(gene_pos_list)):
		while 1:
			pos = int(np.random.random_sample() * L)
			if pos in temp_rep:
				continue
			else:
				temp_rep[pos] = 1
				temp_list.append(pos)
				break
	temp_list = np.array(temp_list)
	temp_list = np.sort(temp_list)
	gene_pos_list = temp_list
	return


def simu_genotype_beta():	# beta (cis-); one more coefficient as the constant item in linear regression
	# to fill in:
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess

	## notes (TODO):
	##	1. an extra item as the constant item, and it's not from Spike and Slab but directly from the Gaussian
	##	2. select some to be non-zero, and all others are zero, with Spike and Slab prior
	##	3. here we don't keep the tissue consistency for these cis- elements, but rather we make them more general; maybe later we can introduce tissue consistency
	##	4. all the parameters in the prior can be changed later on

	global SNP_beta_rep
	global pos_map
	global n_tissue
	global n_gene
	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian

	##==== SNP_beta_rep
	for i in range(n_tissue):
		SNP_beta_rep[i] = {}
		for j in range(n_gene):
			SNP_beta_rep[i][j] = []
			amount = pos_map[j][1] - pos_map[j][0] + 1
			for k in range(amount):
				SNP_beta_rep[i][j].append(0)

			for k in range(len(SNP_beta_rep[i][j])):
				## Spike and Slab prior: first binomial; then 0 or Gaussian			
				# Binomial
				flag = np.random.binomial(1, pi)
				if flag == 1:
					# Gaussian
					beta = np.random.normal(lamb, std)
					SNP_beta_rep[i][j][k] = beta
				else:
					continue

			beta = np.random.normal(lamb, std)
			SNP_beta_rep[i][j].append(beta)
	return


def simu_cell_factor():
	# to fill in:
	##factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
	##beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues

	## notes (TODO):
	##	1. an extra item as the constant item, and it's not from Spike and Slab but directly from the Gaussian
	##	2. select some to be non-zero, and all others are zero, with Spike and Slab prior
	##	3. here we don't keep the tissue consistency for these cis- elements, but rather we make them more general; maybe later we can introduce tissue consistency
	##	4. all the parameters in the prior can be changed later on
	##	5. we set a smaller Slab component, as there are huge amout of trans- SNPs

	global factor_cell_beta_rep
	global n_SNP
	global beta_factor_cell_rep
	global n_factor_cell
	pi = 0.1	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian

	##==== factor_cell_beta_rep
	for i in factor_cell_beta_rep:
		for j in range(n_SNP):
			## Spike and Slab prior: first binomial; then 0 or Gaussian
			# Binomial
			flag = np.random.binomial(1, pi)
			if flag == 1:
				# Gaussian
				beta = np.random.normal(lamb, std)
				factor_cell_beta_rep[i][j] = beta
			else:
				continue
		beta = np.random.normal(lamb, std)
		factor_cell_beta_rep[i][-1] = beta

	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian

	##==== beta_factor_cell_rep
	for i in beta_factor_cell_rep:
		for j in beta_factor_cell_rep[i]:
			for k in range(n_factor_cell):
				## Spike and Slab prior: first binomial; then 0 or Gaussian
				# Binomial
				flag = np.random.binomial(1, pi)
				if flag == 1:
					# Gaussian
					beta = np.random.normal(lamb, std)
					beta_factor_cell_rep[i][j][k] = beta
				else:
					continue
			beta = np.random.normal(lamb, std)
			beta_factor_cell_rep[i][j][-1] = beta
	return



## simulate batch variables for all individuals, and for all samples in different individuals; then concatenate them when using them
def simu_batch_factor():	# variables and beta
	# to fill in:
	##batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
	##batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	##factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
	##beta_factor_batch_rep = {gene:[], ...}

	global n_batch
	global n_batch_individual
	global n_batch_sample
	global n_factor_batch
	global batch_individual_rep
	global batch_tissue_sample_rep
	global factor_batch_beta_rep
	global beta_factor_batch_rep


	##==== batch_individual_rep
	for i in batch_individual_rep:
		for j in range(len(batch_individual_rep[i])):
			value = np.random.random_sample()
			batch_individual_rep[i][j] = value

	##==== batch_tissue_sample_rep
	for i in batch_tissue_sample_rep:
		for j in batch_tissue_sample_rep[i]:
			for k in range(len(batch_tissue_sample_rep[i][j])):	
				value = np.random.random_sample()
				batch_tissue_sample_rep[i][j][k] = value


	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian

	##==== factor_batch_beta_rep
	for i in factor_batch_beta_rep:
		for j in range(n_batch):
			## Spike and Slab prior: first binomial; then 0 or Gaussian
			# Binomial
			flag = np.random.binomial(1, pi)
			if flag == 1:
				# Gaussian
				beta = np.random.normal(lamb, std)
				factor_batch_beta_rep[i][j] = beta
			else:
				continue
		beta = np.random.normal(lamb, std)
		factor_batch_beta_rep[i][-1] = beta


	##==== beta_factor_batch_rep
	for i in beta_factor_batch_rep:
		for j in range(n_factor_batch):
			## Spike and Slab prior: first binomial; then 0 or Gaussian
			# Binomial
			flag = np.random.binomial(1, pi)
			if flag == 1:
				# Gaussian
				beta = np.random.normal(lamb, std)
				beta_factor_batch_rep[i][j] = beta
			else:
				continue
		beta = np.random.normal(lamb, std)
		beta_factor_batch_rep[i][-1] = beta

	return







##=====================
##==== other utilities
## notes: TODO
##	1. xxx
##=====================
def SNP_gene_map():
	# to fill in:
	##pos_map = {}							# {gene:[snp1, snp2], ...}

	global pos_map
	global SNP_pos_list
	global gene_pos_list

	##==== pos_map
	for i in range(len(gene_pos_list)):
		index_gene = 0
		index_start = 0
		index_end = 0
		
		pos_gene = gene_pos_list[i]
		index = 0
		while (SNP_pos_list[index] - pos_gene) < -1000000:
			index += 1
		index_start = index
		index += 1
		while (SNP_pos_list[index] - pos_gene < 1000000):
			index += 1
			if index == len(SNP_pos_list):
				break
		index -= 1
		index_end = index
		pos_map[index_gene] = [index_start, index_end]
	return








if __name__ == '__main__':


	##===========================
	##==== initialize dimensions
	##===========================
	#==== chromosome
	L = 100000000				# 100MB
	#==== individual
	n_individual = 185			# as in GTEx.v.4
	#==== tissue
	n_tissue = 17				# as in GTEx.v.4
	#==== SNP
	n_SNP = 100000				# 100K, so SNP density is 1/1,000bp
	#==== gene
	n_gene = 2000				# 2K, so gene density is 1/50,000bp
	#==== factor_cell
	n_factor_cell = 400			# as evaluated empirically
	#==== factor_batch (analogous to cell factor pathway)
	n_batch = 229
	n_batch_individual = 160		# as in GTEx.v.4
	n_batch_sample = 69			# as in GTEx.v.4
	n_factor_batch = 50			# as evaluated empirically





	##=======================================================
	##==== simulate variables (fill in the dimensions first)
	##=======================================================
	#==== chromosome
	#L = 100000000				# 100MB
	#==== individual
	#n_individual = 185			# as in GTEx.v.4
	#==== tissue
	#n_tissue = 17				# as in GTEx.v.4


	#==== SNP
	#n_SNP = 100000				# 100K, so SNP density is 1/1,000bp
	##SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
	for i in range(n_individual):
		SNP_rep[i] = []
		for j in range(n_SNP):
			SNP_rep[i].append(0)
	##SNP_pos_list = []
	for i in range(n_SNP):
		SNP_pos_list.append(0)
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	# later on, after "SNP gene pos map"
	simu_genotype()


	#==== gene
	#n_gene = 2000				# 2K, so gene density is 1/50,000bp
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals
	for i in range(n_individual):
		gene_rep[i] = {}
		for j in range(n_tissue):
			gene_rep[i][j] = []
			for k in range(n_gene):
				gene_rep[i][j].append(0)
	##gene_pos_list = []
	simu_gene()


	#==== SNP gene pos map
	#pos_map = {}							# {gene:[snp1, snp2], ...}
	SNP_gene_map()
	##SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
	simu_genotype_beta()


	#==== factor_cell
	#n_factor_cell = 400			# as evaluated empirically
	##factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
	for i in range(n_factor_cell):
		factor_cell_beta_rep[i] = []
		for j in range(n_SNP):
			factor_cell_beta_rep[i].append(0)
		factor_cell_beta_rep[i].append(0)
	#beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
	for i in range(n_tissue):
		beta_factor_cell_rep[i] = {}
		for j in range(n_gene):
			beta_factor_cell_rep[i][j] = []
			for k in range(n_factor_cell):
				beta_factor_cell_rep[i][j].append(0)
			beta_factor_cell_rep[i][j].append(0)	# the constant item
	simu_cell_factor()


	#==== factor_batch (analogous to cell factor pathway)
	#n_batch = 229
	#n_batch_individual = 160		# as in GTEx.v.4
	#n_batch_sample = 69			# as in GTEx.v.4
	#n_factor_batch = 400
	##batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
	for i in range(n_individual):
		batch_individual_rep[i] = []
		for j in range(n_batch_individual):
			batch_individual_rep[i].append(0)
	##batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
	for i in range(n_individual):
		batch_tissue_sample_rep[i] = {}
		for j in range(n_tissue):
			batch_tissue_sample_rep[i][j] = []
			for k in range(n_batch_sample):
				batch_tissue_sample_rep[i][j].append(0)
	##factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
	for i in range(n_factor_batch):
		factor_batch_beta_rep[i] = []
		for j in range(n_batch):
			factor_batch_beta_rep[i].append(0)2
		factor_batch_beta_rep[i].append(0)
	##beta_factor_batch_rep = {gene:[], ...}
	for i in range(n_gene):
		beta_factor_batch_rep[i] = []
		for j in range(n_factor_batch):
			beta_factor_batch_rep[i].append(0)
		beta_factor_batch_rep[i].append(0)
	simu_batch_factor()








	##=======================================================
	##==== generate the expression profile through the model
	##=======================================================
	##gene_rep = {individual:{tissue:[], ...}, ...}			# gene expression list for tissue samples for individuals





	











	##======================
	##==== parameter saving
	##======================














