## Input: SNPs (independent sites), 0/1 binary array
## Output: expression level for some genes
## function: generate the coefficients according to the specified graphical model (the neural model in a Bayesian approach), and generate the corresponding expression level for ALL genes
## notes:
##	1. we should simulate tissue specificity, as the modeling takes consideration of that;
##	2. xxx
##	3. xxx


import numpy as np


##==== global variables
## notes: TODO
##	1. shall we simulate Spike and Slab?
##	2. xxx
##=====================
#==== individual
n_individual = 0

#==== SNP
n_SNP = 0
SNP_rep = {}	# {individual:[], xxx:[], ...}			# real value lists, in [0,1], for individuals
SNP_pos_list = []
SNP_beta_rep = {tissue:{gene:[], xxx:[], ...}, xxx:{}, ...}	# cis- SNP beta lists for different genes in tissuess

#==== gene
n_gene = 0
gene_list = []
gene_pos_list = []

#==== SNP gene pos map
pos_map = {}							# {0:[n_1, n_2], 1:[n_1,n_2], ...}

#==== factor_cell
n_factor_cell = 0
factor_cell_beta_rep = {cell factor:[], xxx:[], ...}		# coefficient lists for cell factors
beta_factor_cell_rep = {tissue:{gene:[], xxx:[], ...}, ...}	# coefficient lists of cell factors for genes in tissues

#==== factor_batch (analogous to cell factor pathway)
n_factor_batch = 0						# n_factor_batch = n_factor_batch_individual + n_factor_batch_sample
n_factor_batch_individual = 0
n_factor_batch_sample = 0
factor_batch_beta_rep = {batch factor:[], xxx:[], ...}		# coefficient lists for batch factors

#==== batch variables (will be concatenated from following individual factors and tissue sample factors)
batch_individual_rep = {0:[], 1:[], ...}			# batch variable lists for individuals
batch_tissue_sample_rep = {individual:{tissue:[], ...}, ...}	# batch variable lists for tissue samples in individuals

#==== tissue
n_tissue = 0





##==== simulating variables
## notes: TODO
##	1. simulate variables observed and latent excluding those along the generation process
##	2. xxx
##=====================
def simu_geno():	# variables and beta (cis-)

	return

## simulate batch variables for all individuals, and for all samples in different individuals; then concatenate them when using them
def simu_batch():	# variables and beta

	return

def simu_cell_beta():

	return

def simu_beta_cell_fac():

	return

def simu_beta_batch():

	return








"""
## function: randomly generate the genotypes (minor allele frequency) of n individuals for n_SNP independent sites, value in [0, 0.5, 1]
n = 185
n_SNP = 10000


if __name__ == '__main__':

	file = open("genotype.data", "w")


	for i in range(n):

		array = np.random.randint(3, size=n_SNP)
		array = array * 0.5
		for freq in array:
			file.write(str(freq) + ' ')
		file.write('\n')

	file.close()
"""





if __name__ == '__main__':





	"""
	###==============###
	## generate the positions of associated SNPs with the current gene, fill in index_eQTL
	index_eQTL = np.sort( np.random.permutation(n_SNP)[0: n_eQTL] )  ## sorted index
	power_eQTL = 0.5 * np.random.random_sample((n_eQTL,)) + 0.5



	###==============###
	## generate the TF (cell state function) part --> tissue specificity and consistency
	for i in range(n_Tissue):
		array = np.zeros((n_TF, n_SNP))
		power_TF.append(array)
	for i in range(n_TF):
		## generate the index of associated SNPs, and generate their powers, and sign in matrix across different tissues
		index = np.sort( np.random.permutation(n_SNP)[0: n_TFSNP] )  ## sorted index
		power = 0.5 * np.random.random_sample((n_TFSNP, )) + 0.5

		permutation = np.random.permutation(n_Tissue)

		# keep (n_Tissue - 3) amount of tissues as consistent part
		for j in range(n_Tissue - 3):
			tissue = permutation[j]
			## assign all the power of SNPs into this tissue
			for k in range(len(index)):
				power_TF[tissue][index[k]] = power[k]

		# make (3) tissues show tissue specificity; keep the index as the same, but generate the new power
		for j in range(3):
			power = 0.5 * np.random.random_sample((n_TFSNP, )) + 0.5
			tissue = permutation[n_Tissue - 3 + j]
			## assign all the power of SNPs into this tissue
			for k in range(len(index)):
				power_TF[tissue][index[k]] = power[k]


	## generate the TF (cell state function) part --> the power of TF
	power_afterTF = 0.5 * np.random.random_sample((n_TF, )) + 0.5




	###==============###
	## generate the population part
	"""





	###
	

