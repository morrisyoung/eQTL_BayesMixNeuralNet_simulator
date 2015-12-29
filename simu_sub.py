## functions: simulate given parts in the models
## notes:
##	1. to add the hierarchical regulation (hierarchy code has been finished)
##	2. xxx


import numpy as np
from hierarchy import *




## target:
##	SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
##	SNP_pos_list = []
def simu_genotype(SNP_rep, SNP_pos_list, L):
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
	for i in range(len(SNP_pos_list)):
		SNP_pos_list[i] = temp_list[i]

	return




## target:
##	gene_pos_list
def simu_gene_pos(gene_pos_list, L):
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
	for i in range(len(gene_pos_list)):
		gene_pos_list[i] = temp_list[i]

	return





## target:
##	SNP_beta_rep = {tissue:{gene:[], ...}, ...}			# cis- SNP beta lists for different genes in tissuess
## notes:
##	1. one more coefficient as the constant item in linear regression, which is not from Spike and Slab but directly from the Gaussian
##	2. select some to be non-zero, and all others are zero, with Spike and Slab prior
##	3. here we don't keep the tissue consistency for these cis- elements, but rather we make them more general; try later on
##	4. all the parameters in the prior can be changed later on
def simu_genotype_beta(pos_map, n_tissue, n_gene, SNP_beta_rep):
	# TODO: parameters tunable
	##==============================================================================================================
	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian
	##==============================================================================================================


	# TODO: add the hierarchical regulation
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

			# intercept
			beta = np.random.normal(lamb, std)
			SNP_beta_rep[i][j].append(beta)

	return






## target:
##	factor_cell_beta_rep = {cell factor:[], ...}			# coefficient lists for cell factors
##	beta_factor_cell_rep = {tissue:{gene:[], ...}, ...}		# coefficient lists of cell factors for genes in tissues
## notes:
##	1. an extra item as the constant item, and it's not from Spike and Slab but directly from the Gaussian
##	2. select some to be non-zero, and all others are zero, with Spike and Slab prior
##	3. here we don't keep the tissue consistency for these cis- elements, but rather we make them more general; try later on
##	4. all the parameters in the prior can be tuned later on
##	5. we set a smaller Slab component, as there are handful trans- SNPs with respect to the total candidates
def simu_cell_factor(factor_cell_beta_rep, n_SNP, beta_factor_cell_rep, n_factor_cell):
	# TODO: parameters tunable
	##==============================================================================================================
	pi = 0.1	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian
	##==============================================================================================================

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

		# intercept
		beta = np.random.normal(lamb, std)
		factor_cell_beta_rep[i][-1] = beta

	# TODO: parameters tunable
	##==============================================================================================================
	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian
	##==============================================================================================================


	# TODO: add the hierarchical regulation
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

			# intercept
			beta = np.random.normal(lamb, std)
			beta_factor_cell_rep[i][j][-1] = beta

	return





## target:
##	batch_individual_rep = {individual:[], ...}			# batch variable lists for individuals
##	batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals
##	factor_batch_beta_rep = {batch factor:[], ...}			# coefficient lists for batch factors
##	beta_factor_batch_rep = {gene:[], ...}
## notes:
##	1. this is for: simulate batch variables for all individuals, and for all samples in different individuals; then concatenate them when using them
##	2. simulate variables and the beta
def simu_batch_factor(n_batch, n_batch_individual, n_batch_sample, n_factor_batch, batch_individual_rep, batch_tissue_sample_rep, factor_batch_beta_rep, beta_factor_batch_rep):

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

	# TODO: parameters tunable
	##==============================================================================================================
	pi = 0.5	# for Binomial
	lamb = 0	# for Gaussian
	std = 1		# for Gaussian
	##==============================================================================================================

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

		# intercept
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

		# intercept
		beta = np.random.normal(lamb, std)
		beta_factor_batch_rep[i][-1] = beta

	return



